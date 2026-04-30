using Test
using QAtlas
using QAtlas: TFIM, OBC, Infinite

# Self-validation harness (`test/util/thermodynamic_identities.jl`)
# applied to TFIM at every BC where the required quantities are
# implemented.  Two identities ship in DEFAULT_IDENTITIES:
#
#   - Gibbs:                 ε = f + T·s
#   - SpecificHeat-from-ε:   c_v = -β² ∂ε/∂β   (ForwardDiff)
#
# Both use `Energy(:per_site)` (granted by the granularity dispatch
# from PR #124) so per-site/total mix-ups in TFIM would surface here
# as :fail rather than a silent wrong number.

@testset "TFIM ε = f + T·s and c_v = -β² ∂ε/∂β  — OBC(N=8)" begin
    model = TFIM(; J=1.0, h=0.5)
    βs = [0.5, 1.0, 2.0]
    results = verify_thermodynamic_identities(model, OBC(8); βs=βs)

    # 2 identities × 3 βs
    @test length(results) == 6
    @test all(r.status === :pass for r in results)

    # Spot-check that both identities actually ran (didn't all skip)
    @test any(occursin("Gibbs", r.identity) for r in results)
    @test any(occursin("c_v", r.identity) for r in results)

    # Numerical tightness — TFIM has closed-form thermodynamics, so the
    # residuals should sit well below the harness's default `(1e-8, 1e-10)`.
    for r in results
        @test r.abs_err < 1e-8
    end
end

@testset "TFIM ε = f + T·s and c_v = -β² ∂ε/∂β  — Infinite()" begin
    model = TFIM(; J=1.0, h=0.5)
    βs = [0.5, 1.0, 2.0]
    results = verify_thermodynamic_identities(model, Infinite(); βs=βs)

    @test length(results) == 6
    @test all(r.status === :pass for r in results)

    # The Infinite() Energy + thermal observables go through QuadGK; tighter
    # tolerance than OBC dense ED but still within the harness threshold.
    for r in results
        @test r.abs_err < 1e-7
    end
end

@testset "TFIM at the critical point h = J — identities still hold" begin
    # The critical point is the most numerically delicate slice (gap closes,
    # entropy peaks).  Self-validation must still hold there.
    model = TFIM(; J=1.0, h=1.0)
    βs = [1.0, 5.0]  # avoid extreme β to keep dense-ED + QuadGK stable
    for bc in (OBC(8), Infinite())
        results = verify_thermodynamic_identities(model, bc; βs=βs)
        @test all(r.status === :pass for r in results)
    end
end

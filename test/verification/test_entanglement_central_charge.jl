# ─────────────────────────────────────────────────────────────────────────────
# Verification: central charge from entanglement entropy
#
# For a 1D critical system with OBC, the Calabrese-Cardy formula gives:
#
#   S(l) = (c/6) ln[(2N/π) sin(πl/N)] + s₁
#
# where c is the central charge and s₁ is non-universal.
#
# Note on the Peschel correlation-matrix method: the TFIM with σ^z σ^z
# convention maps to free fermions only after a Kramers-Wannier duality.
# The BdG eigenvectors give the DUAL fermion correlators, not the spin-
# basis correlators. Correct spin-basis Peschel requires handling the
# duality boundary conditions, which is deferred to future work.
# For now, we use full ED (exact for small N, O(2^N) cost).
#
# References:
#   P. Calabrese, J. Cardy, J. Stat. Mech. 0406, P06002 (2004).
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Lattice2D, LinearAlgebra, Test

include("../util/spinhalf_ed.jl")

# Extract c from S(l) profile using OBC Calabrese-Cardy (c/6)
function extract_central_charge_obc(Ss::Vector{Float64}, ls, N::Int)
    ξs = [log((2N / π) * sin(π * l / N)) for l in ls]
    n = length(ls)
    ξ̄ = sum(ξs) / n
    S̄ = sum(Ss) / n
    slope = sum((ξs[i] - ξ̄) * (Ss[i] - S̄) for i in 1:n) /
            sum((ξs[i] - ξ̄)^2 for i in 1:n)
    return 6 * slope  # OBC: S = (c/6) ξ + const
end

@testset "Entanglement entropy → central charge (ED)" begin

    @testset "TFIM at h = J → c = 1/2" begin
        c_exact = Float64(QAtlas.fetch(Universality(:Ising), CriticalExponents(); d=2).c)
        J = 1.0

        # Multiple N for finite-size scaling of the extraction
        cs_extracted = Dict{Int,Float64}()
        for N in [10, 12, 14]
            lat = build_lattice(Square, N, 1; boundary=OpenAxis())
            H = build_tfim(lat, J, J)
            ψ0 = eigen(Symmetric(H)).vectors[:, 1]

            # Skip 2 boundary sites on each side (lattice artifacts)
            ls = collect(3:(N - 3))
            Ss = [entanglement_entropy(ψ0, l, N) for l in ls]
            cs_extracted[N] = extract_central_charge_obc(Ss, ls, N)
        end

        # Larger N gives better c (finite-size corrections decrease)
        @test abs(cs_extracted[14] - c_exact) < abs(cs_extracted[10] - c_exact)

        # N=14 should be within 10% of c = 1/2
        @test cs_extracted[14] ≈ c_exact rtol = 0.10

        # S(l) structural checks at N=12
        N = 12
        lat = build_lattice(Square, N, 1; boundary=OpenAxis())
        H = build_tfim(lat, J, J)
        ψ0 = eigen(Symmetric(H)).vectors[:, 1]
        Ss = [entanglement_entropy(ψ0, l, N) for l in 1:(N - 1)]

        # Symmetric: S(l) ≈ S(N-l)
        for l in 2:(N ÷ 2 - 1)
            @test Ss[l] ≈ Ss[N - l] rtol = 0.01
        end

        # Maximal near center
        l_max = argmax(Ss)
        @test abs(l_max - N / 2) <= 1
    end

    @testset "TFIM: area law away from criticality" begin
        N = 10
        lat = build_lattice(Square, N, 1; boundary=OpenAxis())

        # Disordered phase (h ≫ J): product state → S ≈ 0
        H_dis = build_tfim(lat, 1.0, 10.0)
        ψ_dis = eigen(Symmetric(H_dis)).vectors[:, 1]
        S_dis = entanglement_entropy(ψ_dis, N ÷ 2, N)

        # Critical: S is larger
        H_crit = build_tfim(lat, 1.0, 1.0)
        ψ_crit = eigen(Symmetric(H_crit)).vectors[:, 1]
        S_crit = entanglement_entropy(ψ_crit, N ÷ 2, N)

        @test S_dis < 0.1
        @test S_crit > S_dis * 5
    end

    @testset "Heisenberg chain → c = 1" begin
        J = 1.0
        N = 12
        lat = build_lattice(Square, N, 1; boundary=OpenAxis())
        H = build_spinhalf_heisenberg(lat, J)
        ψ0 = eigen(Symmetric(H)).vectors[:, 1]

        # Even l only (suppress SU(2) alternating correction)
        ls_even = collect(4:2:(N - 4))
        Ss = [entanglement_entropy(ψ0, l, N) for l in ls_even]
        c_heis = extract_central_charge_obc(Ss, ls_even, N)

        @test c_heis ≈ 1.0 rtol = 0.20

        # c(Heisenberg) > c(TFIM)
        H_tfim = build_tfim(lat, J, J)
        ψ_tfim = eigen(Symmetric(H_tfim)).vectors[:, 1]
        ls_all = collect(3:(N - 3))
        Ss_tfim = [entanglement_entropy(ψ_tfim, l, N) for l in ls_all]
        c_tfim = extract_central_charge_obc(Ss_tfim, ls_all, N)

        @test c_heis > c_tfim
    end
end

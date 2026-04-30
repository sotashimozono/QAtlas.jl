using Test
using QAtlas
using QAtlas:
    Energy,
    OBC,
    PBC,
    Infinite,
    TFIM,
    XXZ1D,
    Heisenberg1D,
    S1Heisenberg1D,
    KitaevHoneycomb,
    native_energy_granularity,
    fetch

# Tests for the parameterized `Energy{G}` dispatch surface introduced
# alongside the `native_energy_granularity` trait.  Goal: pin down
# (i) backward-compat of bare `Energy()`, (ii) explicit-granularity
# round-trips through the generic conversion fallback, and (iii) the
# Kitaev (2D, per-site native) corner where the `:total` fallback is
# expected to error.

@testset "Energy{G} constructors and aliases" begin
    @test Energy() isa Energy{:natural}
    @test Energy(:natural) isa Energy{:natural}
    @test Energy(:total) isa Energy{:total}
    @test Energy(:per_site) isa Energy{:per_site}
    @test_throws ErrorException Energy(:bogus)
    @test_throws ErrorException Energy{:bogus}()
end

@testset "native_energy_granularity declarations" begin
    @test native_energy_granularity(TFIM(; J=1.0, h=1.0), OBC(8)) === :total
    @test native_energy_granularity(TFIM(; J=1.0, h=1.0), Infinite()) === :per_site
    @test native_energy_granularity(XXZ1D(; J=1.0, Δ=0.0), OBC(8)) === :total
    @test native_energy_granularity(XXZ1D(; J=1.0, Δ=0.0), Infinite()) === :per_site
    @test native_energy_granularity(Heisenberg1D(), OBC(8)) === :total
    @test native_energy_granularity(S1Heisenberg1D(; J=1.0), OBC(4)) === :total
    @test native_energy_granularity(KitaevHoneycomb(), OBC(0)) === :per_site
    @test native_energy_granularity(KitaevHoneycomb(), PBC(0)) === :per_site
    @test native_energy_granularity(KitaevHoneycomb(), Infinite()) === :per_site
end

@testset "Energy() :natural router routes to native granularity" begin
    # TFIM OBC: :natural → :total
    let m = TFIM(; J=1.0, h=0.5), N = 8, β = 1.0
        @test fetch(m, Energy(), OBC(N); beta=β) == fetch(m, Energy(:total), OBC(N); beta=β)
    end
    # TFIM Infinite: :natural → :per_site
    let m = TFIM(; J=1.0, h=0.5), β = 1.0
        @test fetch(m, Energy(), Infinite(); beta=β) ==
            fetch(m, Energy(:per_site), Infinite(); beta=β)
    end
    # Kitaev (per-site at every BC)
    let m = KitaevHoneycomb()
        @test fetch(m, Energy(), Infinite()) == fetch(m, Energy(:per_site), Infinite())
        @test fetch(m, Energy(), OBC(0); Lx=3, Ly=3) ==
            fetch(m, Energy(:per_site), OBC(0); Lx=3, Ly=3)
        @test fetch(m, Energy(), PBC(0); Lx=3, Ly=3) ==
            fetch(m, Energy(:per_site), PBC(0); Lx=3, Ly=3)
    end
end

@testset "Energy(:per_site) at finite BC = Energy(:total) / N (1D conversion)" begin
    # TFIM
    for N in (4, 6, 8), β in (0.5, 1.0, 2.0)
        m = TFIM(; J=1.0, h=0.7)
        e_tot = fetch(m, Energy(:total), OBC(N); beta=β)
        e_ps = fetch(m, Energy(:per_site), OBC(N); beta=β)
        @test e_ps ≈ e_tot / N rtol=1e-12
    end
    # XXZ1D
    for N in (4, 6), β in (0.5, 1.0), Δ in (-1.0, 0.0, 1.0)
        m = XXZ1D(; J=1.0, Δ=Δ)
        e_tot = fetch(m, Energy(:total), OBC(N); beta=β)
        e_ps = fetch(m, Energy(:per_site), OBC(N); beta=β)
        @test e_ps ≈ e_tot / N rtol=1e-12
    end
    # Heisenberg1D (delegates to XXZ1D)
    for N in (4, 6), β in (0.5, 1.0)
        e_tot = fetch(Heisenberg1D(), Energy(:total), OBC(N); beta=β, J=1.5)
        e_ps = fetch(Heisenberg1D(), Energy(:per_site), OBC(N); beta=β, J=1.5)
        @test e_ps ≈ e_tot / N rtol=1e-12
    end
    # S1Heisenberg1D
    for N in (3, 4), β in (0.5, 1.0)
        m = S1Heisenberg1D(; J=1.0)
        e_tot = fetch(m, Energy(:total), OBC(N); beta=β)
        e_ps = fetch(m, Energy(:per_site), OBC(N); beta=β)
        @test e_ps ≈ e_tot / N rtol=1e-12
    end
end

@testset "Energy(:total) at Infinite is undefined" begin
    # Generic conversion only covers OBC/PBC; at Infinite no model
    # registers Energy{:total} natively (total energy diverges in the
    # thermodynamic limit), so the request falls through to QAtlas's
    # informative catch-all `fetch` (`src/core/type.jl`), which raises
    # an `ErrorException`.
    @test_throws ErrorException fetch(
        TFIM(; J=1.0, h=1.0), Energy(:total), Infinite(); beta=1.0
    )
    @test_throws ErrorException fetch(XXZ1D(; J=1.0, Δ=0.0), Energy(:total), Infinite())
    @test_throws ErrorException fetch(KitaevHoneycomb(), Energy(:total), Infinite())
end

@testset "Energy(:total) on 2D Kitaev: generic fallback fails (expected)" begin
    # Kitaev declares native :per_site at OBC/PBC.  The generic :total
    # fallback would call `_bc_size(bc, kwargs)`, but Kitaev encodes the
    # system size in `Lx, Ly` kwargs rather than `bc.N`, so the conversion
    # is intentionally unsupported and surfaces the standard "N
    # unspecified" error.  This test pins the behaviour so a future Kitaev
    # PR that adds 2D system_size support is forced to update both the
    # implementation and this expectation.
    @test_throws ErrorException fetch(KitaevHoneycomb(), Energy(:total), OBC(0); Lx=3, Ly=3)
    @test_throws ErrorException fetch(KitaevHoneycomb(), Energy(:total), PBC(0); Lx=3, Ly=3)
end

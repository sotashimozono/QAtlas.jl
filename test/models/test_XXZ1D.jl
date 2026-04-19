using QAtlas, Test

@testset "XXZ1D — dispatch & construction" begin
    m = XXZ1D()
    @test m isa QAtlas.AbstractQAtlasModel
    @test m.J == 1.0
    @test m.Δ == 0.0

    m2 = XXZ1D(; J=2.5, Δ=0.7)
    @test m2.J == 2.5
    @test m2.Δ == 0.7

    # Symbol alias resolves to :XXZ1D canonical
    @test QAtlas.canonicalize_model(Val(:XXZ)) === :XXZ1D
    @test QAtlas.canonicalize_model(Val(:xxz)) === :XXZ1D
    @test QAtlas.canonicalize_model(Val(:xxz1d)) === :XXZ1D
end

@testset "XXZ1D — known closed-form points (J=1, Infinite)" begin
    # Δ = 0 (XX, free fermion): e₀/J = -1/π
    e_xx = QAtlas.fetch(XXZ1D(; J=1.0, Δ=0.0), Energy(), Infinite())
    @test e_xx ≈ -1 / π atol = 1e-10

    # Δ = 1 (AF Heisenberg, Hulthén 1938): e₀/J = 1/4 - ln 2
    e_af = QAtlas.fetch(XXZ1D(; J=1.0, Δ=1.0), Energy(), Infinite())
    @test e_af ≈ 0.25 - log(2.0) atol = 1e-10

    # Δ = -1 (FM): e₀/J = -1/4
    e_fm = QAtlas.fetch(XXZ1D(; J=1.0, Δ=-1.0), Energy(), Infinite())
    @test e_fm ≈ -0.25 atol = 1e-14

    # GroundStateEnergyDensity alias returns the same value.
    e_gs = QAtlas.fetch(XXZ1D(; J=1.0, Δ=1.0), GroundStateEnergyDensity(), Infinite())
    @test e_gs ≈ e_af
end

@testset "XXZ1D — Energy outside the three canonical Δ is NaN + warning" begin
    # The general-Δ Bethe-ansatz integral is deferred; for now XXZ1D
    # advertises only Δ ∈ {-1, 0, 1} and returns NaN + a warning elsewhere.
    for Δ in (-0.5, 0.3, 0.7, 1.5, 2.0)
        e = @test_logs (:warn, r"general-Δ") QAtlas.fetch(
            XXZ1D(; J=1.0, Δ=Δ), Energy(), Infinite()
        )
        @test isnan(e)
    end
end

@testset "XXZ1D — Luttinger parameter K" begin
    @test QAtlas.fetch(XXZ1D(; J=1.0, Δ=0.0), LuttingerParameter(), Infinite()) ≈ 1.0 atol =
        1e-12
    @test QAtlas.fetch(XXZ1D(; J=1.0, Δ=1.0), LuttingerParameter(), Infinite()) ≈ 0.5 atol =
        1e-12

    # Monotone decreasing from Δ = -1 (K → ∞) to Δ = 1 (K = 1/2)
    ks = [
        QAtlas.fetch(XXZ1D(; J=1, Δ=Δ), LuttingerParameter(), Infinite()) for
        Δ in -0.9:0.2:0.9
    ]
    @test all(diff(ks) .< 0)
end

@testset "XXZ1D — LuttingerVelocity (+ SpinWaveVelocity alias)" begin
    # Canonical values at Δ=0 (XX) and Δ=1 (AF Heisenberg)
    @test QAtlas.fetch(XXZ1D(; J=1.0, Δ=0.0), LuttingerVelocity(), Infinite()) ≈ 1.0 atol =
        1e-12
    @test QAtlas.fetch(XXZ1D(; J=1.0, Δ=1.0), LuttingerVelocity(), Infinite()) ≈ π / 2 atol =
        1e-12

    # J scaling linear
    @test QAtlas.fetch(XXZ1D(; J=3.0, Δ=0.0), LuttingerVelocity(), Infinite()) ≈ 3.0

    # SpinWaveVelocity() === LuttingerVelocity at the type level
    @test SpinWaveVelocity === LuttingerVelocity
    @test SpinWaveVelocity() isa LuttingerVelocity

    u1 = QAtlas.fetch(XXZ1D(; J=1.0, Δ=0.5), LuttingerVelocity(), Infinite())
    u2 = QAtlas.fetch(XXZ1D(; J=1.0, Δ=0.5), SpinWaveVelocity(), Infinite())
    @test u1 === u2
end

@testset "XXZ1D — central charge (critical regime only)" begin
    for Δ in (-0.9, -0.5, 0.0, 0.5, 0.99)
        @test QAtlas.fetch(XXZ1D(; J=1.0, Δ=Δ), CentralCharge(), Infinite()) == 1.0
    end
    # Outside the critical regime we return NaN + a warning.
    @test isnan(
        @test_logs (:warn, r"critical regime") QAtlas.fetch(
            XXZ1D(; J=1.0, Δ=1.5), CentralCharge(), Infinite()
        )
    )
end

@testset "XXZ1D — gapped regime is NaN (general-Δ deferred)" begin
    # |Δ| > 1: gapped; NaN + warning.  This is the same branch as the
    # generic-Δ path above (single deferred-implementation warning);
    # the test spells it out for documentation value.
    e = @test_logs (:warn, r"general-Δ") QAtlas.fetch(
        XXZ1D(; J=1.0, Δ=2.0), Energy(), Infinite()
    )
    @test isnan(e)
end

@testset "XXZ1D — legacy Symbol dispatch routes to new API" begin
    # Deprecation log on every new (model, quantity) pair
    e_sym = QAtlas.fetch(:XXZ, :energy, Infinite(); J=1.0, Δ=0.0)
    @test e_sym ≈ -1 / π atol = 1e-10

    u_sym = QAtlas.fetch(:XXZ, :spin_wave_velocity, Infinite(); J=1.0, Δ=1.0)
    @test u_sym ≈ π / 2 atol = 1e-12

    u_fv = QAtlas.fetch(:XXZ, :fermi_velocity, Infinite(); J=1.0, Δ=0.0)
    @test u_fv ≈ 1.0 atol = 1e-12

    K_sym = QAtlas.fetch(:XXZ, :luttinger_parameter, Infinite(); J=1.0, Δ=0.0)
    @test K_sym ≈ 1.0 atol = 1e-12
end

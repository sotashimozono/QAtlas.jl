# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: Onsager critical temperature + Yang magnetization
#
# Lattice-independent checks of the closed-form IsingSquare results.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "IsingSquare — Onsager T_c" begin
    Tc = QAtlas.fetch(IsingSquare(), CriticalTemperature())
    @test Tc ≈ 2 / log(1 + sqrt(2)) atol = 1e-14
    @test Tc ≈ 2.269185314213022 atol = 1e-10  # known numerical value

    # J scaling: T_c(J) = J · T_c(1)
    for J in (0.5, 2.0, 3.7)
        @test QAtlas.fetch(IsingSquare(), CriticalTemperature(); J=J) ≈ J * Tc rtol = 1e-14
    end

    # Relation to sinh: at T_c, sinh(2J/T_c) = 1
    β_c = 1 / Tc
    @test sinh(2 * β_c) ≈ 1.0 atol = 1e-14
end

@testset "IsingSquare — Yang spontaneous magnetization" begin
    Tc = QAtlas.fetch(IsingSquare(), CriticalTemperature())
    β_c = 1 / Tc

    @testset "Special values" begin
        # T = 0 (β → ∞): M = 1 (fully ordered)
        @test QAtlas.fetch(IsingSquare(), SpontaneousMagnetization(); β=1000.0) ≈ 1.0 atol =
            1e-10

        # T = T_c: M = 0
        @test QAtlas.fetch(IsingSquare(), SpontaneousMagnetization(); β=β_c) == 0.0

        # T > T_c: M = 0
        @test QAtlas.fetch(IsingSquare(), SpontaneousMagnetization(); β=β_c * 0.5) == 0.0

        # T slightly below T_c: M is positive (but NOT small — the exponent
        # β = 1/8 is so small that M ~ (1 - T/T_c)^{1/8} ≈ 0.01^{1/8} ≈ 0.7
        # even at 1% below T_c). Just verify positivity and < 1.
        M_near = QAtlas.fetch(IsingSquare(), SpontaneousMagnetization(); β=β_c * 1.01)
        @test 0 < M_near < 1
    end

    @testset "M is monotonically increasing as T → 0 (β increasing)" begin
        betas = [β_c * 1.1, β_c * 1.5, β_c * 2, β_c * 5, β_c * 20]
        Ms = [QAtlas.fetch(IsingSquare(), SpontaneousMagnetization(); β=β) for β in betas]
        for k in 1:(length(Ms) - 1)
            @test Ms[k] < Ms[k + 1]
        end
    end

    @testset "Critical exponent β = 1/8 near T_c" begin
        # M ~ (T_c - T)^{1/8} as T → T_c⁻
        # log M / log(T_c - T) → 1/8
        # Test at two temperatures near T_c
        δT1 = 0.001
        δT2 = 0.01
        T1 = Tc - δT1
        T2 = Tc - δT2
        M1 = QAtlas.fetch(IsingSquare(), SpontaneousMagnetization(); β=1 / T1)
        M2 = QAtlas.fetch(IsingSquare(), SpontaneousMagnetization(); β=1 / T2)
        # Effective exponent: log(M1/M2) / log(δT1/δT2) ≈ β = 1/8
        β_eff = log(M1 / M2) / log(δT1 / δT2)
        @test β_eff ≈ 1 / 8 rtol = 0.05  # within 5% of exact
    end

    @testset "J scaling" begin
        β = 0.5
        M1 = QAtlas.fetch(IsingSquare(), SpontaneousMagnetization(); β=β, J=1.0)
        # M(β, J) = M(βJ=const), so M(β=0.5, J=2) = M(β=1, J=1)
        M2 = QAtlas.fetch(IsingSquare(), SpontaneousMagnetization(); β=1.0, J=1.0)
        M_scaled = QAtlas.fetch(IsingSquare(), SpontaneousMagnetization(); β=0.5, J=2.0)
        @test M_scaled ≈ M2 atol = 1e-14
    end
end

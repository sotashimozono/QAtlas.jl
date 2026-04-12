# ─────────────────────────────────────────────────────────────────────────────
# Standalone test: universality class critical exponents
#
# Verify exact (rational) exponents for each universality class and
# check standard scaling relations between them.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Test

@testset "Universality: 2D Ising critical exponents" begin
    e = QAtlas.fetch(Ising2D(), CriticalExponents())

    # Exact rational values
    @test e.β == 1 // 8
    @test e.ν == 1 // 1
    @test e.γ == 7 // 4
    @test e.η == 1 // 4
    @test e.δ == 15 // 1
    @test e.α == 0 // 1
    @test e.c == 1 // 2

    # Standard scaling relations (must hold exactly for rational exponents)
    @test e.α + 2 * e.β + e.γ == 2               # Rushbrooke
    @test e.γ == e.β * (e.δ - 1)                  # Widom
    @test e.γ == e.ν * (2 - e.η)                  # Fisher
    @test 2 - e.α == 2 * e.ν                      # Josephson (d=2: dν = 2-α)
end

@testset "Universality: KPZ 1+1D scaling exponents" begin
    e = QAtlas.fetch(KPZ1D(), CriticalExponents())

    @test e.β_growth == 1 // 3
    @test e.α_rough == 1 // 2
    @test e.z == 3 // 2

    # KPZ scaling relations
    @test e.α_rough + e.z == 2                    # Galilean invariance
    @test e.β_growth == e.α_rough / e.z           # β = α / z
end

@testset "Universality: Mean-field critical exponents" begin
    e = QAtlas.fetch(MeanField(), CriticalExponents())

    @test e.β == 1 // 2
    @test e.ν == 1 // 2
    @test e.γ == 1 // 1
    @test e.η == 0 // 1
    @test e.δ == 3 // 1
    @test e.α == 0 // 1

    # Scaling relations (hold at mean-field level)
    @test e.α + 2 * e.β + e.γ == 2               # Rushbrooke
    @test e.γ == e.β * (e.δ - 1)                  # Widom
    @test e.γ == e.ν * (2 - e.η)                  # Fisher
end


@testset "E8 Spectrum Logic" begin
    # 1. 基本的な構造のテスト
    @testset "Structure" begin
        masses = QAtlas.fetch(E8(), E8Spectrum(), Infinite())
        @test length(masses) == 8
        @test masses[1] == 1.0
        @test all(masses .> 0)
        @test issorted(masses) # 質量は昇順であるべき
    end

    # 2. m₂/m₁ = golden ratio
    @testset "Golden Ratio Relationship" begin
        masses = QAtlas.fetch(E8(), E8Spectrum(), Infinite())
        ϕ = (1 + sqrt(5)) / 2
        @test masses[2] ≈ ϕ atol = 1e-15
    end

    # 3. alias — deliberately reaches through the legacy Symbol shim to
    # verify that `:mass_ratios`, `:E8_masses`, `:mass_ratio` all
    # canonicalise to `:E8_spectrum` and forward to `E8Spectrum`.  The
    # shim emits a one-shot `@info` per (model, quantity) pair; wrap
    # every legacy call with `@test_logs` so the deprecation notices do
    # not leak into CI output.
    @testset "Aliases" begin
        expected = QAtlas.fetch(E8(), E8Spectrum(), Infinite())

        @test @test_logs (:info, r"symbol-dispatch") QAtlas.fetch(:E8, :mass_ratios) ==
            expected
        @test @test_logs (:info, r"symbol-dispatch") QAtlas.fetch(:E8, :E8_masses) ==
            expected
        @test @test_logs (:info, r"symbol-dispatch") QAtlas.fetch(:E8, :mass_ratio) ==
            expected
    end

    @testset "Type Stability" begin
        @test @inferred(QAtlas.fetch(Model(:E8), Quantity(:E8_spectrum), Infinite())) isa
            Vector{Float64}
    end
end

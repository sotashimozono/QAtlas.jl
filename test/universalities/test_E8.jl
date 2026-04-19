@testset "E8 Spectrum Logic" begin
    # ───────────────────────────────────────────────────────────────────
    # 1. Structural sanity: 8 masses, ordered, m₁ = 1, all positive.
    # ───────────────────────────────────────────────────────────────────
    @testset "Structure" begin
        masses = QAtlas.fetch(E8(), E8Spectrum(), Infinite())
        @test length(masses) == 8
        @test masses[1] == 1.0
        @test all(masses .> 0)
        @test issorted(masses)  # 質量は昇順であるべき
    end

    # ───────────────────────────────────────────────────────────────────
    # 2. Golden-ratio identity m₂/m₁ = φ — hallmark of E₈ symmetry,
    #    measured experimentally in Coldea et al., Science 327, 177 (2010).
    # ───────────────────────────────────────────────────────────────────
    @testset "Golden Ratio Relationship" begin
        masses = QAtlas.fetch(E8(), E8Spectrum(), Infinite())
        ϕ = (1 + sqrt(5)) / 2
        @test masses[2] ≈ ϕ atol = 1e-15
    end

    # ───────────────────────────────────────────────────────────────────
    # 3. Literature cross-check: every stored mass must equal the
    #    closed-form Zamolodchikov (1989) / Delfino (2004) expression
    #    to machine precision. The expressions below are derived
    #    independently in `docs/src/calc/e8-mass-spectrum-derivation.md`
    #    (Perron-Frobenius eigenvector of the E₈ Cartan matrix + Dorey
    #    fusion rules). Any sign flip, coefficient drift, or angle
    #    typo in `get_e8_mass_ratios()` is caught here.
    #
    #    Historical note: v0.13.4 and earlier silently returned
    #    `m₇, m₈` with a factor-of-2 error (the `4 · ϕ² · cos(…)`
    #    spelling of the original implementation collapses to
    #    `2 · ϕ² · cos(…)` once the definition of φ is expanded).
    #    This testset would have failed on that code — kept as a
    #    regression guard.
    # ───────────────────────────────────────────────────────────────────
    @testset "Closed-form expressions (Zamolodchikov 1989, Delfino 2004)" begin
        masses = QAtlas.fetch(E8(), E8Spectrum(), Infinite())
        ϕ = 2 * cos(π / 5)

        @test masses[1] == 1.0
        @test masses[2] ≈ ϕ atol = 1e-14
        @test masses[3] ≈ 2 * cos(π / 30) atol = 1e-14
        @test masses[4] ≈ 2 * cos(7π / 30) * ϕ atol = 1e-14
        @test masses[5] ≈ 2 * cos(2π / 15) * ϕ atol = 1e-14
        @test masses[6] ≈ 2 * cos(π / 30) * ϕ atol = 1e-14
        @test masses[7] ≈ 2 * cos(7π / 30) * ϕ^2 atol = 1e-14
        @test masses[8] ≈ 2 * cos(2π / 15) * ϕ^2 atol = 1e-14
    end

    # ───────────────────────────────────────────────────────────────────
    # 4. Decimal-level cross-check against the tabulated literature
    #    values (Delfino 2004 review, eq. (4.14); Zamolodchikov 1989
    #    Table 2).  These are the six-digit numbers every subsequent
    #    paper and textbook quotes.  A trivial coefficient bug in the
    #    mass formula would shift them by O(1), so a tight atol is
    #    effective.
    # ───────────────────────────────────────────────────────────────────
    @testset "Tabulated decimal values (Delfino 2004, eq. 4.14)" begin
        masses = QAtlas.fetch(E8(), E8Spectrum(), Infinite())
        expected = [
            1.000000,
            1.618034,
            1.989044,
            2.404867,
            2.956295,
            3.218340,
            3.891157,
            4.783386,
        ]
        for i in 1:8
            @test masses[i] ≈ expected[i] atol = 1e-5
        end
    end

    # ───────────────────────────────────────────────────────────────────
    # 5. Fusion-rule cross-check: several pairs of masses are linked by
    #    the on-shell bootstrap fusion `a × b → c`, which forces
    #    multiplicative identities between mass ratios. Two clean ones
    #    (Zamolodchikov 1989, §5; Delfino 2004 Table 2):
    #
    #      m₆ / (m₂ m₃)  = 1         (fusion 2 × 3 → 6)
    #      m₇ / (m₂ m₄)  = 1         (fusion 2 × 4 → 7)
    #      m₈ / (m₂ m₅)  = 1         (fusion 2 × 5 → 8)
    #
    #    These have to hold exactly (up to floating-point round-off) in
    #    the closed-form spectrum, independent of the derivation path,
    #    so they are a strong internal cross-check that the stored
    #    values are consistent with the integrable-field-theory
    #    bootstrap.
    # ───────────────────────────────────────────────────────────────────
    @testset "Fusion-rule multiplicative identities" begin
        m = QAtlas.fetch(E8(), E8Spectrum(), Infinite())
        @test m[6] ≈ m[2] * m[3] atol = 1e-14
        @test m[7] ≈ m[2] * m[4] atol = 1e-14
        @test m[8] ≈ m[2] * m[5] atol = 1e-14
    end

    # ───────────────────────────────────────────────────────────────────
    # 6. Legacy Symbol-dispatch routes the aliases `:mass_ratios`,
    #    `:E8_masses`, `:mass_ratio` all canonicalise to `:E8_spectrum`
    #    and forward to `E8Spectrum`.  The shim emits a one-shot
    #    `@info` per (model, quantity) pair; wrap every legacy call
    #    with `@test_logs` so the deprecation notices do not leak into
    #    CI output.
    # ───────────────────────────────────────────────────────────────────
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

using QAtlas, Lattice2D, LinearAlgebra, Test

# The full-ED reference needs the spin-ED utilities loaded by runtests.jl;
# when this file is included from the test harness they are already in scope.

@testset "TFIM VonNeumannEntropy — Peschel vs full ED at small N" begin
    # N = 10 → 2^N = 1024 full state vector; the SVD reference is
    # cheap and the Peschel covariance is a 20 × 20 matrix.
    N = 10
    for (J, h) in ((1.0, 0.3), (1.0, 1.0), (1.0, 3.0))
        lat = build_lattice(Square, N, 1; boundary=OpenAxis())
        H = build_tfim(lat, J, h)
        ψ0 = eigen(Symmetric(H)).vectors[:, 1]
        model = TFIM(; J=J, h=h)
        for ℓ in 1:(N - 1)
            S_ed = entanglement_entropy(ψ0, ℓ, N)
            S_peschel = QAtlas.fetch(model, VonNeumannEntropy(), OBC(N); ℓ=ℓ)
            @test S_peschel ≈ S_ed atol = 1e-8
        end
    end
end

@testset "TFIM VonNeumannEntropy — symmetry and structure" begin
    # S(ℓ) = S(N − ℓ) for pure states (ground state); equivalent to
    # saying the Peschel covariance of the two halves gives the same
    # spectrum.
    N = 16
    model = TFIM(; J=1.0, h=1.0)  # critical, where entanglement is richest
    Ss = [QAtlas.fetch(model, VonNeumannEntropy(), OBC(N); ℓ=ℓ) for ℓ in 1:(N - 1)]
    for ℓ in 1:(N ÷ 2 - 1)
        @test Ss[ℓ] ≈ Ss[N - ℓ] rtol = 1e-10
    end
    @test argmax(Ss) == N ÷ 2 || argmax(Ss) == N ÷ 2 - 1

    # Disordered phase (h ≫ J): S ≪ 1 for any ℓ (area law; here exactly 0
    # in the h → ∞ product state, near-zero at h = 10 J).
    model_dis = TFIM(; J=1.0, h=10.0)
    S_dis = QAtlas.fetch(model_dis, VonNeumannEntropy(), OBC(N); ℓ=N ÷ 2)
    @test S_dis < 0.1

    # Ordered phase far from criticality: also area-law (ignoring the
    # near-degenerate Z_2 tunnel splitting at this finite N).
    model_ord = TFIM(; J=1.0, h=0.1)
    S_ord = QAtlas.fetch(model_ord, VonNeumannEntropy(), OBC(N); ℓ=N ÷ 2)
    @test S_ord < 1.1  # upper-bounded by ln 2 (Z_2 cat state) + finite-size noise
end

@testset "TFIM VonNeumannEntropy — Calabrese–Cardy log scaling at h = J" begin
    # S(ℓ, N) = (c/6) log( (2N/π) sin(πℓ/N) ) + s₁ for an OBC chain
    # at the Ising CFT critical point with c = 1/2.  At N = 100 the
    # Peschel result should fit this form with c extracted within 5 %.
    N = 100
    model = TFIM(; J=1.0, h=1.0)
    ℓs = collect(10:10:(N - 10))
    Ss = [QAtlas.fetch(model, VonNeumannEntropy(), OBC(N); ℓ=ℓ) for ℓ in ℓs]
    ξs = [log((2N / π) * sin(π * ℓ / N)) for ℓ in ℓs]
    n = length(ℓs)
    ξ̄ = sum(ξs) / n
    S̄ = sum(Ss) / n
    slope = sum((ξs[i] - ξ̄) * (Ss[i] - S̄) for i in 1:n) /
            sum((ξs[i] - ξ̄)^2 for i in 1:n)
    c_est = 6 * slope
    @test c_est ≈ 0.5 rtol = 0.05
end

@testset "TFIM VonNeumannEntropy — thermal beta argument" begin
    # At β = Inf the thermal covariance reduces to the ground-state
    # covariance; the two call paths must agree to machine precision.
    N = 20
    model = TFIM(; J=1.0, h=1.0)
    ℓ = N ÷ 2
    S_gs = QAtlas.fetch(model, VonNeumannEntropy(), OBC(N); ℓ=ℓ)
    S_inf = QAtlas.fetch(model, VonNeumannEntropy(), OBC(N); ℓ=ℓ, beta=Inf)
    @test S_gs ≈ S_inf atol = 1e-12

    # Finite β < ∞: thermal entropy + entanglement entropy is larger
    # than the T = 0 value (thermal mixing adds classical entropy).
    S_thermal = QAtlas.fetch(model, VonNeumannEntropy(), OBC(N); ℓ=ℓ, beta=1.0)
    @test S_thermal > S_gs
end

@testset "TFIM VonNeumannEntropy — legacy Symbol dispatch" begin
    N = 10
    S_new = QAtlas.fetch(TFIM(; J=1.0, h=1.0), VonNeumannEntropy(), OBC(N); ℓ=N ÷ 2)
    S_legacy = QAtlas.fetch(
        :TFIM, :entanglement_entropy, OBC(); N=N, J=1.0, h=1.0, ℓ=N ÷ 2
    )
    @test S_legacy ≈ S_new atol = 1e-12

    # Symbol aliases :ee, :EE, :S_vN, :EntanglementEntropy all canonicalise
    # to :entanglement_entropy.
    for q in (:ee, :EE, :S_vN, :EntanglementEntropy)
        @test QAtlas.fetch(:TFIM, q, OBC(); N=N, J=1.0, h=1.0, ℓ=N ÷ 2) ≈ S_new atol = 1e-12
    end
end

@testset "TFIM VonNeumannEntropy — input validation" begin
    @test_throws ArgumentError QAtlas.fetch(
        TFIM(; J=1.0, h=1.0), VonNeumannEntropy(), OBC(10); ℓ=0
    )
    @test_throws ArgumentError QAtlas.fetch(
        TFIM(; J=1.0, h=1.0), VonNeumannEntropy(), OBC(10); ℓ=10
    )
end

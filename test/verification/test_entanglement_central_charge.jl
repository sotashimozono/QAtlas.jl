# ─────────────────────────────────────────────────────────────────────────────
# Verification: central charge extraction from entanglement entropy
#
# For a 1D critical system with OBC and N sites, the Calabrese-Cardy
# formula (2004) gives the bipartite entanglement entropy:
#
#   S(l) = (c/3) ln[(2N/π) sin(πl/N)] + s₁
#
# where c is the central charge and s₁ is a non-universal constant.
# By computing S(l) for several subsystem sizes l and fitting against
# the conformal coordinate ξ(l) = ln[(2N/π) sin(πl/N)], we extract c.
#
# Cross-check:
#   Source A: Universality(:Ising) → c = 1/2  (CFT, BPZ 1984)
#   Source B: TFIM ED ground state → S(l) → c  (Calabrese-Cardy 2004)
#
# This is the DEFINITIVE cross-verification of the central charge,
# directly from the ground-state wavefunction.
# ─────────────────────────────────────────────────────────────────────────────

using QAtlas, Lattice2D, LinearAlgebra, Test

include("../util/spinhalf_ed.jl")

"""
    extract_central_charge(ψ, N; bc=:OBC) -> Float64

Extract the central charge c from the entanglement entropy profile
S(l) of a ground state using the Calabrese-Cardy formula.

For **OBC** (one entanglement cut):
    S(l) = (c/6) ln[(2N/π) sin(πl/N)] + s₁    →  c = 6 × slope

For **PBC** (two entanglement cuts):
    S(l) = (c/3) ln[(N/π) sin(πl/N)] + s₁     →  c = 3 × slope

Uses linear regression of S(l) vs the conformal coordinate ξ(l).

# References
    P. Calabrese, J. Cardy, J. Stat. Mech. 0406, P06002 (2004), Eqs. (7) and (19).
"""
function extract_central_charge(ψ::AbstractVector, N::Int; bc::Symbol=:OBC)
    ls = 1:(N - 1)
    Ss = [entanglement_entropy(ψ, l, N) for l in ls]

    if bc == :OBC
        ξs = [log((2N / π) * sin(π * l / N)) for l in ls]
        prefactor = 6  # S = (c/6) ξ + const for OBC
    else  # PBC
        ξs = [log((N / π) * sin(π * l / N)) for l in ls]
        prefactor = 3  # S = (c/3) ξ + const for PBC
    end

    n = length(ls)
    ξ̄ = sum(ξs) / n
    S̄ = sum(Ss) / n
    slope = sum((ξs[i] - ξ̄) * (Ss[i] - S̄) for i in 1:n) / sum((ξs[i] - ξ̄)^2 for i in 1:n)
    return prefactor * slope
end

@testset "Entanglement entropy → central charge" begin
    @testset "TFIM at criticality → c = 1/2" begin
        c_exact = Float64(QAtlas.fetch(Universality(:Ising), CriticalExponents(); d=2).c)
        J = 1.0
        h = J  # critical point

        # Use N = 12 OBC chain (2^12 = 4096 dim Hilbert space)
        N = 12
        lat = build_lattice(Square, N, 1; boundary=OpenAxis())
        H = build_tfim(lat, J, h)
        F = eigen(Symmetric(H))
        ψ0 = F.vectors[:, 1]  # ground state

        c_extracted = extract_central_charge(ψ0, N; bc=:OBC)
        # With the correct OBC Calabrese-Cardy prefactor (c/6 instead
        # of c/3 for PBC), N=12 should give c ≈ 0.5 within ~10%.
        @test c_extracted ≈ c_exact rtol = 0.15

        # S(l) should be maximal near the center l ≈ N/2
        Ss = [entanglement_entropy(ψ0, l, N) for l in 1:(N - 1)]
        l_max = argmax(Ss)
        @test abs(l_max - N / 2) <= 1  # peak near center

        # S should be symmetric: S(l) ≈ S(N-l)
        for l in 1:(N ÷ 2 - 1)
            @test Ss[l] ≈ Ss[N - l] rtol = 0.05
        end
    end

    @testset "TFIM away from criticality → area law (low S)" begin
        J = 1.0
        N = 10
        lat = build_lattice(Square, N, 1; boundary=OpenAxis())

        # Deep in disordered phase: h ≫ J → product state → S ≈ 0
        H_dis = build_tfim(lat, J, 10.0)
        ψ_dis = eigen(Symmetric(H_dis)).vectors[:, 1]
        S_center = entanglement_entropy(ψ_dis, N ÷ 2, N)
        @test S_center < 0.1  # nearly zero (area law)

        # Compare to critical: S should be much larger
        H_crit = build_tfim(lat, J, J)
        ψ_crit = eigen(Symmetric(H_crit)).vectors[:, 1]
        S_crit = entanglement_entropy(ψ_crit, N ÷ 2, N)
        @test S_crit > S_center * 5  # critical ≫ gapped
    end

    @testset "Heisenberg chain → c = 1 (Luttinger liquid)" begin
        # The 1D Heisenberg AFM chain is a c = 1 CFT (free boson /
        # Luttinger liquid). This is double the TFIM value.
        #
        # CAVEAT: The Heisenberg chain has well-known ALTERNATING
        # corrections (−1)^l f(l) to the Calabrese-Cardy formula, arising
        # from the SU(2) symmetry. These make naive all-l regression less
        # accurate. We use only EVEN l values to suppress the oscillation.
        J = 1.0
        N = 12
        lat = build_lattice(Square, N, 1; boundary=OpenAxis())
        H = build_spinhalf_heisenberg(lat, J)
        F = eigen(Symmetric(H))
        ψ0 = F.vectors[:, 1]

        # Use only even l values (suppress alternating correction from SU(2))
        ls_even = 2:2:(N - 2)
        Ss = [entanglement_entropy(ψ0, l, N) for l in ls_even]
        ξs = [log((2N / π) * sin(π * l / N)) for l in ls_even]
        n = length(ls_even)
        ξ̄ = sum(ξs) / n
        S̄ = sum(Ss) / n
        slope =
            sum((ξs[i] - ξ̄) * (Ss[i] - S̄) for i in 1:n) / sum((ξs[i] - ξ̄)^2 for i in 1:n)
        c_heis = 6 * slope  # OBC: c = 6 × slope

        # c = 1 for Heisenberg (OBC + finite-size + alternating corrections
        # → 20% tolerance; even-l filtering suppresses oscillation)
        @test c_heis ≈ 1.0 rtol = 0.20

        # c(Heisenberg) should be larger than c(TFIM)
        lat_tfim = build_lattice(Square, N, 1; boundary=OpenAxis())
        H_tfim = build_tfim(lat_tfim, J, J)
        ψ_tfim = eigen(Symmetric(H_tfim)).vectors[:, 1]
        c_tfim = extract_central_charge(ψ_tfim, N; bc=:OBC)

        @test c_heis > c_tfim  # c=1 > c=1/2
    end
end

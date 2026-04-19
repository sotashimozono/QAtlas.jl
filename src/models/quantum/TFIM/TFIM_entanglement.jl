# ─────────────────────────────────────────────────────────────────────────────
# Transverse Field Ising Model — ground / thermal entanglement entropy via
# Peschel's correlation-matrix method.
#
# For a contiguous block A = {site 1, …, site ℓ} of the N-site OBC TFIM, the
# reduced density matrix ρ_A of a Gaussian BdG state (ground state or thermal
# state) is itself Gaussian.  Let
#
#   Σ_A = Σ[1:2ℓ, 1:2ℓ]                   (2ℓ × 2ℓ real antisymmetric)
#
# be the Majorana covariance (see `_majorana_thermal_covariance` in
# `TFIM_dynamics.jl`) restricted to the 2ℓ Majoranas on sites 1..ℓ.  The
# 2ℓ eigenvalues of the Hermitian matrix `i Σ_A` come in ± pairs
# {±ν_1, …, ±ν_ℓ} with 0 ≤ ν_k ≤ 1; equivalently ν_k are the singular
# values of Σ_A.  Each ν_k is the covariance eigenvalue of an
# independent Bogoliubov mode of ρ_A, with occupation
# n_k = (1 ± ν_k)/2 and per-mode von Neumann entropy
#
#   s(ν_k) = -[(1-ν_k)/2 · ln((1-ν_k)/2) + (1+ν_k)/2 · ln((1+ν_k)/2)].
#
# Summing,
#
#   S_A = Σ_{k=1}^{ℓ} s(ν_k).                                         (1)
#
# JW-convention remark.  The existing `_majorana_ham` uses the σˣ-string
# JW convention (see the header of `TFIM_dynamics.jl`): the Majorana
# pair (γ_{2i-1}, γ_{2i}) is local to spin site i.  For a *contiguous*
# bipartition A = {sites 1..ℓ} the JW unitary factorises across the A/B
# boundary up to a parity factor on A that commutes with ρ_A, so the
# fermion entropy S_A^{(f)} equals the spin entropy S_A^{(s)}.  (Peschel
# 2003; Fagotti–Calabrese 2010 for the factorisation argument.)
#
# Cost.  O(ℓ³) from the Hermitian eigendecomposition of `i Σ_A`.  For
# N = 200, ℓ = 100 the calculation takes ~10 ms.  The full-ED baseline
# (SVD of the 2^N ground-state vector, `test_entanglement_central_charge.jl`)
# is O(2^N) memory and O(4^N) time and breaks down at N ≳ 16.
#
# References:
# - I. Peschel, J. Phys. A 36, L205 (2003), eq. (9).
# - G. Vidal, J. I. Latorre, E. Rico, A. Kitaev, Phys. Rev. Lett. 90,
#   227902 (2003), §III.
# - J. I. Latorre, E. Rico, G. Vidal, Quantum Inf. Comput. 4, 48 (2004),
#   §2.
# - M. Fagotti, P. Calabrese, Phys. Rev. Lett. 104, 227203 (2010) —
#   JW-factorisation argument for contiguous blocks.
# ─────────────────────────────────────────────────────────────────────────────

using LinearAlgebra: eigvals, Hermitian

"""
    _peschel_mode_entropy(ν::Float64) -> Float64

Per-mode von Neumann entropy

    s(ν) = -[(1-ν)/2 · ln((1-ν)/2) + (1+ν)/2 · ln((1+ν)/2)]

as a function of the Majorana covariance singular value `ν ∈ [0, 1]`.
Numerically stable at `ν → 0` (`s → ln 2`) and `ν → 1` (`s → 0`).
"""
function _peschel_mode_entropy(ν::Float64)::Float64
    ν = clamp(abs(ν), 0.0, 1.0)
    p⁻ = (1 - ν) / 2
    p⁺ = (1 + ν) / 2
    s⁻ = p⁻ > 1e-15 ? -p⁻ * log(p⁻) : 0.0
    s⁺ = p⁺ > 1e-15 ? -p⁺ * log(p⁺) : 0.0
    return s⁻ + s⁺
end

"""
    fetch(model::TFIM, ::VonNeumannEntropy, bc::OBC;
          ℓ::Int, beta::Float64 = Inf, kwargs...) -> Float64

Von Neumann entanglement entropy of the first `ℓ` spins of the N-site
OBC TFIM in the thermal state at inverse temperature `beta` (or the
ground state when `beta = Inf`), computed by Peschel's correlation-
matrix method — see equation (1) in the file header.

- `N = _bc_size(bc, kwargs)` (read from `OBC(N)` or legacy `kwargs[:N]`).
- `ℓ` must satisfy `1 ≤ ℓ ≤ N - 1`.
- Cost is `O(ℓ³)`; for typical `N = 200, ℓ = 100` this runs in a few
  milliseconds, whereas the full-ED SVD baseline scales as `O(4^N)`.

The result matches the full-ED reference at every small `N` (verified
to 1e-10 in `test/models/test_TFIM_entanglement.jl`).

See full derivation: [JW reduction of the TFIM](../../../../docs/src/calc/jw-tfim-bdg.md).
"""
function fetch(
    model::TFIM, ::VonNeumannEntropy, bc::OBC; ℓ::Int, beta::Float64=Inf, kwargs...
)
    N = _bc_size(bc, kwargs)
    1 ≤ ℓ ≤ N - 1 || throw(
        ArgumentError(
            "VonNeumannEntropy: ℓ must satisfy 1 ≤ ℓ ≤ N - 1; got ℓ = $ℓ, N = $N."
        ),
    )
    hmat = _majorana_ham(N, model.J, model.h)
    Σ = _majorana_thermal_covariance(hmat, beta)
    Σ_A = Σ[1:(2ℓ), 1:(2ℓ)]
    # Eigenvalues of `i Σ_A` are real (Hermitian) and come in ± ν pairs.
    # Ascending sort places the ℓ non-negative ν's in the upper half.
    λ = eigvals(Hermitian(im .* Σ_A))
    S = 0.0
    @inbounds for k in (ℓ + 1):(2ℓ)
        S += _peschel_mode_entropy(λ[k])
    end
    return S
end

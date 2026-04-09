# ─────────────────────────────────────────────────────────────────────────────
# Transverse Field Ising Model — exact solutions
#
# Hamiltonian:
#   H = -J Σᵢ σᶻᵢσᶻᵢ₊₁  -  h Σᵢ σˣᵢ
#
# Solved exactly via Jordan-Wigner + Bogoliubov-de Gennes (BdG) transformation.
# The quadratic fermion Hamiltonian has quasiparticle energies Λₙ > 0, giving:
#
#   ⟨H⟩(β) = -Σₙ (Λₙ/2) tanh(β Λₙ / 2)
#
# Aliases: :TFIM, :TransverseFieldIsingModel
# ─────────────────────────────────────────────────────────────────────────────

using LinearAlgebra: eigvals, Symmetric
using QuadGK: quadgk

# ═══════════════════════════════════════════════════════════════════════════════
# Internal: BdG quasiparticle spectrum (OBC, finite N)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _tfim_bdg_spectrum(N, J, h) -> Vector{Float64}

Return the N positive BdG quasiparticle energies Λₙ > 0 for the OBC TFIM
with N sites, Ising coupling J, and transverse field h.

The 2N×2N BdG matrix is:
    H_BdG = [[A, B]; [-B, -A]]
where A (tridiagonal, symmetric) encodes hopping + onsite energy,
and B (antisymmetric) encodes the pairing terms from JW transformation.

    A_{ii}   = 2h
    A_{i,i±1} = -J
    B_{i,i+1} = +J,  B_{i+1,i} = -J
"""
function _tfim_bdg_spectrum(N::Int, J::Float64, h::Float64)::Vector{Float64}
    A = zeros(N, N)
    for i in 1:N
        A[i, i] = 2h
    end
    for i in 1:(N - 1)
        A[i, i + 1] = -J
        A[i + 1, i] = -J
    end

    B = zeros(N, N)
    for i in 1:(N - 1)
        B[i, i + 1] = J
        B[i + 1, i] = -J
    end

    H_bdg = [A B; -B -A]
    vals = eigvals(Symmetric(H_bdg))
    return sort!(filter(v -> v > 1e-10, vals))
end

# ═══════════════════════════════════════════════════════════════════════════════
# Energy: OBC finite-N
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(model::Model{:TFIM}, ::Quantity{:energy}, ::OBC; beta, betas) -> Float64 or Vector{Float64}

Total energy ⟨H⟩(β) for the OBC TFIM with N sites.

Required model params: `N` (Int), `J` (Float64), `h` (Float64, transverse field)

- `beta::Float64`: return scalar ⟨H⟩(β)
- `betas::AbstractVector{Float64}`: return vector, reusing spectrum (O(N³) once)
- no keyword: return ground-state energy E₀ = -Σₙ Λₙ/2  (β → ∞)

Uses the exact BdG formula:  ⟨H⟩ = -Σₙ (Λₙ/2) tanh(β Λₙ / 2)
"""
function fetch(
    model::Model{:TFIM},
    ::Quantity{:energy},
    ::OBC;
    beta::Union{Float64,Nothing}=nothing,
    betas::Union{AbstractVector{Float64},Nothing}=nothing,
)
    N = Int(model.params[:N])
    J = Float64(model.params[:J])
    h = Float64(model.params[:h])
    Λ = _tfim_bdg_spectrum(N, J, h)
    if betas !== nothing
        return Float64[-sum(λ -> (λ / 2) * tanh(β * λ / 2), Λ) for β in betas]
    elseif beta !== nothing
        return -sum(λ -> (λ / 2) * tanh(beta * λ / 2), Λ)
    else
        # Ground state: β → ∞, tanh(β Λ/2) → 1
        return -sum(Λ) / 2
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# Energy: thermodynamic limit (PBC / Infinite)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(model::Model{:TFIM}, ::Quantity{:energy}, ::Infinite; beta, betas) -> Float64 or Vector{Float64}

Energy *per site* ⟨H⟩/N in the thermodynamic limit (PBC, N → ∞).

    ε(β) = -(1/π) ∫₀^π dk  Λ(k)/2 · tanh(β Λ(k) / 2)

where the PBC dispersion is  Λ(k) = 2√(J² + h² - 2Jh cos k).

- `beta::Float64`: return scalar ε(β)
- `betas::AbstractVector{Float64}`: return vector
- no keyword: return ground-state energy per site (β → ∞)

Uses adaptive Gauss-Kronrod quadrature (QuadGK).
"""
function fetch(
    model::Model{:TFIM},
    ::Quantity{:energy},
    ::Infinite;
    beta::Union{Float64,Nothing}=nothing,
    betas::Union{AbstractVector{Float64},Nothing}=nothing,
)
    J = Float64(model.params[:J])
    h = Float64(model.params[:h])
    _energy_at_beta =
        β -> begin
            result, _ = quadgk(
                k -> begin
                    Λk = 2sqrt(J^2 + h^2 - 2J * h * cos(k))
                    (Λk / 2) * tanh(β * Λk / 2)
                end, 0.0, π; rtol=1e-10
            )
            -(1 / π) * result
        end
    if betas !== nothing
        return [_energy_at_beta(β) for β in betas]
    elseif beta !== nothing
        return _energy_at_beta(beta)
    else
        # Ground state: β → ∞
        return _energy_at_beta(1e6)
    end
end

# ═══════════════════════════════════════════════════════════════════════════════
# Central charge (critical point h = J)
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(model::Model{:TFIM}, ::Quantity{:central_charge}, ::Infinite) -> Float64

Central charge of the TFIM critical point (h = J): c = 1/2 (Ising CFT).
Returns NaN if not at the critical point (|h/J - 1| > 1e-6).
"""
function fetch(model::Model{:TFIM}, ::Quantity{:central_charge}, ::Infinite)
    J = Float64(model.params[:J])
    h = Float64(model.params[:h])
    if abs(h / J - 1.0) > 1e-6
        @warn "TFIM central charge c=1/2 only at the critical point h=J. Got h/J=$(h/J)."
        return NaN
    end
    return 0.5
end

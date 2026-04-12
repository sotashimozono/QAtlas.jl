# ─────────────────────────────────────────────────────────────────────────────
# Heisenberg — spin-1/2 antiferromagnetic Heisenberg model
#
# Hamiltonian:
#   H = J Σ_{⟨i,j⟩} S_i · S_j
#
# where S_i are spin-1/2 operators and ⟨i,j⟩ runs over nearest-neighbor
# pairs. For J > 0 the ground state is a singlet.
#
# ─────────────────────────────────────────────────────────────────────────────
# Dimer (N = 2) — total-spin analysis
#
#   Two spin-1/2 degrees of freedom combine into S_tot = 0 ⊕ S_tot = 1.
#   The dimer Hamiltonian is
#
#     H = J S_1 · S_2 = (J/2)·[(S_1 + S_2)² − S_1² − S_2²]
#       = (J/2)·[S_tot² − 3/2]
#
#   giving eigenvalues
#
#     S_tot = 0 (singlet):   E_s = −3J/4
#     S_tot = 1 (triplet):   E_t = +J/4  (three-fold degenerate)
#
#   Singlet–triplet gap: Δ = E_t − E_s = J.
# ─────────────────────────────────────────────────────────────────────────────

# ═══════════════════════════════════════════════════════════════════════════════
# Dispatch tags
# ═══════════════════════════════════════════════════════════════════════════════

"""
    Heisenberg1D

Dispatch tag for the spin-1/2 antiferromagnetic Heisenberg model on a 1D
chain (or more generally any finite spin-1/2 cluster). Hamiltonian:

    H = J Σ_{⟨i,j⟩} S_i · S_j,   spin-1/2, J > 0 antiferromagnetic
"""
struct Heisenberg1D end

"""
    ExactSpectrum

Dispatch tag for the full sorted eigenvalue spectrum of a finite model.
"""
struct ExactSpectrum end

# ═══════════════════════════════════════════════════════════════════════════════
# fetch: exact spectrum for small N
# ═══════════════════════════════════════════════════════════════════════════════

"""
    fetch(::Heisenberg1D, ::ExactSpectrum; N, J=1.0) -> Vector{Float64}

Return the sorted exact spectrum of the spin-1/2 Heisenberg Hamiltonian
on an N-site cluster with a single nearest-neighbor bond per pair.

Currently only `N = 2` (the dimer) is implemented, with the
closed-form spectrum

    [-3J/4, J/4, J/4, J/4]

corresponding to one singlet (S_tot = 0) and a three-fold degenerate
triplet (S_tot = 1).

Longer chains (4-site PBC with E₀ = −2J, etc.) are tracked as future
additions in `dev/1_physical-verification/todo.md`.

# Arguments
- `N::Int`: number of spin-1/2 sites (only `N = 2` supported for now)
- `J::Real`: Heisenberg coupling constant (default 1.0; J > 0 AFM)

# References
    For the dimer derivation see any standard textbook, e.g.
    A. Auerbach, "Interacting Electrons and Quantum Magnetism" (1994), §2.
"""
function fetch(::Heisenberg1D, ::ExactSpectrum; N::Int, J::Real=1.0)
    if N == 2
        return sort([-3J / 4, J / 4, J / 4, J / 4])
    end
    return error(
        "Heisenberg1D exact spectrum: only N=2 (dimer) implemented; got N=$N. " *
        "Longer chains are tracked in dev/1_physical-verification/todo.md.",
    )
end

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
on an N-site chain or ring with boundary condition `bc`.

# Supported cases

- **N=2, bc=:OBC** (dimer): `[-3J/4, J/4, J/4, J/4]`
  (singlet E_s = -3J/4, triplet E_t = J/4, three-fold degenerate).
- **N=4, bc=:PBC** (4-site ring): `[-2J, -J×3, 0×7, +J×5]`
  Ground state E₀ = -2J (unique singlet). The ferromagnetic quintet
  sits at E = +J. The full degeneracy structure is:
  1 singlet + 1 triplet + (1 singlet + 2 triplets at E=0) + 1 quintet.

# Arguments
- `N::Int`: number of spin-1/2 sites
- `J::Real`: Heisenberg coupling constant (default 1.0; J > 0 AFM)
- `bc::Symbol`: boundary condition, `:OBC` (default) or `:PBC`

# References
    A. Auerbach, "Interacting Electrons and Quantum Magnetism" (1994), §2.
    H. Bethe, Z. Physik 71, 205 (1931).
"""
function fetch(::Heisenberg1D, ::ExactSpectrum; N::Int, J::Real=1.0, bc::Symbol=:OBC)
    if N == 2 && bc == :OBC
        return sort([-3J / 4, J / 4, J / 4, J / 4])
    elseif N == 4 && bc == :PBC
        # Exact spectrum of the 4-site PBC Heisenberg ring H = J Σ S_i·S_{i+1}:
        #   E = -2J ×1 (singlet), -J ×3 (triplet), 0 ×7, +J ×5 (quintet)
        return sort([
            -2J,
            -J,
            -J,
            -J,
            zero(J),
            zero(J),
            zero(J),
            zero(J),
            zero(J),
            zero(J),
            zero(J),
            J,
            J,
            J,
            J,
            J,
        ])
    end
    return error(
        "Heisenberg1D exact spectrum: only (N=2, OBC) and (N=4, PBC) " *
        "are implemented; got (N=$N, bc=$bc).",
    )
end

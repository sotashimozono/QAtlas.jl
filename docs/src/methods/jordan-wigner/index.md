# Jordan-Wigner Transformation

## What It Is

The Jordan-Wigner (JW) transformation is an exact mapping between
spin-1/2 operators on a 1D chain and spinless fermion operators. It
converts certain spin Hamiltonians into quadratic (free) fermion
problems that can be solved exactly.

## Definition

$$c_i = \left(\prod_{j<i}\sigma^z_j\right)\sigma^-_i, \qquad c^\dagger_i = \left(\prod_{j<i}\sigma^z_j\right)\sigma^+_i$$

The product $\prod_{j<i}\sigma^z_j$ (the **JW string**) ensures
fermionic anti-commutation: $\{c_i, c^\dagger_j\} = \delta_{ij}$.

## When It Applies

| Condition | Satisfied? | Consequence |
| --------- | ---------- | ----------- |
| 1D geometry | Required | JW string is well-defined only in 1D |
| Nearest-neighbor interactions | Helpful | JW string cancels for NN terms |
| No $\sigma^x\sigma^x$ with $\sigma^z\sigma^z$ mixing | Helpful | Avoids quartic fermion terms |

## Limitations

- **Strictly 1D**: In 2D+, the JW string creates non-local
  interactions. Extensions exist (2D JW, Kitaev) but are more complex.
- **Convention dependence**: $\sigma^z\sigma^z$ and $\sigma^x\sigma^x$
  conventions produce different fermion Hamiltonians — one may be
  free while the other is interacting. See
  [TFIM calculation](../../calc/jw-tfim-bdg.md) for how the
  Kramers-Wannier duality resolves this.

## Applications in QAtlas

| Model | Calculation note | Result |
| ----- | ---------------- | ------ |
| [TFIM](../../models/quantum/tfim.md) | [JW-TFIM-BdG](../../calc/jw-tfim-bdg.md) | BdG quasiparticle spectrum |
| Heisenberg | (interacting after JW — not directly solvable) | Bethe ansatz instead |
| Kitaev chain | (planned) | Topological phase boundary |

## References

- P. Jordan, E. Wigner, Z. Physik **47**, 631 (1928).
- E. Lieb, T. Schultz, D. Mattis, Ann. Phys. **16**, 407 (1961).

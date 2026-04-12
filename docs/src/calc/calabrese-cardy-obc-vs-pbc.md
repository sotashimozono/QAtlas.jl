# Calabrese-Cardy: OBC vs PBC Prefactor

## Setup

A 1+1-dimensional critical system (CFT with central charge $c$) on a
lattice of $N$ sites with UV cutoff $a$. The entanglement entropy of a
contiguous subsystem of $l$ sites is computed by the replica trick
and a conformal mapping.

## Calculation

### PBC (periodic boundary conditions)

The system lives on a ring of circumference $N$. A bipartition into
$[1, l]$ and $[l+1, N]$ produces **two** entanglement cuts (the
subsystem has two boundary points with the complement).

The conformal mapping takes the plane to a cylinder of circumference
$N$. The two-point function of twist operators at separation $l$
on the cylinder gives (Calabrese-Cardy 2004, Eq. (7)):

$$S_{\mathrm{PBC}} = \frac{c}{3}\ln\!\left[\frac{N}{\pi a}\sin\!\left(\frac{\pi l}{N}\right)\right] + s_1$$

The prefactor $c/3$ reflects the two cuts, each contributing $c/6$.

### OBC (open boundary conditions)

The system lives on a segment $[1, N]$. A bipartition at site $l$
produces **one** entanglement cut (only one boundary point between
subsystem and complement, since the chain ends are physical
boundaries).

The conformal mapping takes the upper half-plane to a strip of width
$N$. The one-point function of the twist operator near the boundary
gives (Calabrese-Cardy 2004, Eq. (19)):

$$S_{\mathrm{OBC}} = \frac{c}{6}\ln\!\left[\frac{2N}{\pi a}\sin\!\left(\frac{\pi l}{N}\right)\right] + s_1'$$

The prefactor $c/6$ reflects the single cut. The factor $2N$ (instead
of $N$) arises from the method of images in the strip geometry.

### Physical reason for the factor-of-2 difference

|                    | PBC (cylinder)       | OBC (strip)            |
| ------------------ | -------------------- | ---------------------- |
| Entanglement cuts  | 2                    | 1                      |
| Twist operators    | 2-point function     | 1-point function       |
| Conformal geometry | plane $\to$ cylinder | half-plane $\to$ strip |
| Prefactor          | $c/3$                | $c/6$                  |

Each entanglement cut contributes $c/6$ to the logarithmic
coefficient. PBC has two cuts, OBC has one.

### Non-universal constants

The additive constants $s_1$ (PBC) and $s_1'$ (OBC) are
non-universal and depend on the microscopic details of the model.
They differ between PBC and OBC even for the same model.

## Result

$$
\boxed{S_{\mathrm{PBC}} = \frac{c}{3}\ln\!\left[\frac{N}{\pi a}\sin\frac{\pi l}{N}\right] + s_1,
  \qquad
  S_{\mathrm{OBC}} = \frac{c}{6}\ln\!\left[\frac{2N}{\pi a}\sin\frac{\pi l}{N}\right] + s_1'}
$$

The ratio of logarithmic prefactors is exactly 2, corresponding to
the number of entanglement cuts.

## References

- P. Calabrese, J. Cardy, J. Stat. Mech. (2004) P06002; Eq. (7) for PBC, Eq. (19) for OBC.
- C. Holzhey, F. Larsen, F. Wilczek, Nucl. Phys. B **424**, 443 (1994) — earlier derivation for PBC.

## Used by

- [Calabrese-Cardy Method](../methods/calabrese-cardy/index.md)
- [Entanglement Entropy Verification](../verification/entanglement.md)

# E8 Mass Spectrum: Derivation from Bootstrap Equations

## Setup

The Ising field theory perturbed by the magnetic operator $\sigma$
(see [magnetic perturbation](ising-cft-magnetic-perturbation.md)) is an
integrable massive field theory with **eight stable particles**. Their
mass ratios are universal constants determined by the $E_8$ Lie algebra.

**Goal**: derive the mass ratios $m_n / m_1$ for $n = 1, \ldots, 8$.

## Step 1: Integrability and Elastic S-matrix

Zamolodchikov (1989) showed that the perturbed theory has eight
conserved charges of spins $s = 1, 7, 11, 13, 17, 19, 23, 29$
(the exponents of $E_8$). The existence of these charges implies:

1. The number of particles is conserved in every scattering event
   (no particle production or annihilation).
2. The set of momenta is conserved (purely elastic scattering).
3. The $n$-body S-matrix factorizes into products of two-body S-matrices
   (Yang-Baxter integrability).

## Step 2: Bootstrap Equations

The two-body S-matrix $S_{ab}(\theta)$ (where $\theta$ is the rapidity
difference) must satisfy:

1. **Unitarity**: $S_{ab}(\theta) S_{ab}(-\theta) = 1$
2. **Crossing symmetry**: $S_{ab}(\theta) = S_{a\bar{b}}(i\pi - \theta)$
3. **Bootstrap (fusion)**: if particle $c$ appears as a bound state of
   $a$ and $b$, then
   $$m_c^2 = m_a^2 + m_b^2 + 2m_a m_b \cos u_{ab}^c$$
   where $u_{ab}^c$ is the fusion angle, and the S-matrix satisfies the
   bootstrap equation
   $$S_{dc}(\theta) = S_{da}(\theta + i\bar{u}_{bc}^a) S_{db}(\theta - i\bar{u}_{ac}^b)$$

## Step 3: E8 Root System and Mass Ratios

The solution to the bootstrap equations corresponds to the $E_8$
affine Toda field theory. The eight particle masses are determined by
the **Perron-Frobenius eigenvector** of the $E_8$ Cartan matrix.

Equivalently, the mass ratios can be expressed in terms of the golden
ratio $\varphi = 2\cos(\pi/5) = (1+\sqrt{5})/2$:

## Result

| Particle | Mass ratio $m_n / m_1$ | Exact expression | Numerical value |
| -------- | ---------------------- | ---------------- | --------------- |
| $m_1$ | $1$ | $1$ | $1.000$ |
| $m_2$ | $\varphi$ | $2\cos(\pi/5)$ | $1.618$ |
| $m_3$ | | $2\cos(\pi/30)$ | $1.989$ |
| $m_4$ | | $2\varphi\cos(7\pi/30)$ | $2.405$ |
| $m_5$ | | $2\varphi\cos(2\pi/15)$ | $2.956$ |
| $m_6$ | | $2\varphi\cos(\pi/30)$ | $3.218$ |
| $m_7$ | | $4\varphi^2\cos(7\pi/30)$ | $3.891$ |
| $m_8$ | | $4\varphi^2\cos(2\pi/15)$ | $4.783$ |

### Key features

- $m_2/m_1 = \varphi$ (the **golden ratio**) — the most famous prediction.
- Only the first three particles ($m_1, m_2, m_3$) are below the
  two-particle threshold $2m_1$, so they are absolutely stable.
  Particles $m_4$ through $m_8$ are above threshold but are stabilized
  by integrability (elastic scattering prevents decay).
- The ratios involve only trigonometric functions of rational multiples
  of $\pi$, reflecting the root system of $E_8$.

### Fusion structure

The lightest particles satisfy the fusion rules:

$$m_1 + m_1 \to m_2, \qquad m_1 + m_2 \to m_3, \qquad m_1 + m_3 \to m_4$$

The fusion angles are:

$$m_2^2 = 2m_1^2(1 + \cos(\pi/5)), \quad \text{giving } m_2 = 2m_1\cos(\pi/10)$$

Wait — let me be more precise. The actual relation is:

$$m_2 = 2m_1\cos(\pi/5) = \varphi \cdot m_1$$

which follows from the fusing $1 + 1 \to 2$ with fusion angle
$u_{11}^2 = \pi/5$.

## Physical Significance

The E8 mass spectrum is one of the most striking examples of
**emergent mathematical structure** in physics:

- The $E_8$ Lie algebra is the largest exceptional simple Lie algebra
  (rank 8, dimension 248).
- It appears here not as a symmetry of the Hamiltonian (the symmetry
  is only $\mathbb{Z}_2$) but as a structure of the **S-matrix**.
- The appearance of the golden ratio $\varphi$ in a quantum spin chain
  is a non-trivial prediction that has been experimentally verified.

## Experimental Confirmation

Coldea et al. (2010) studied the quasi-1D Ising ferromagnet CoNb₂O₆
near its quantum critical point ($B_c \approx 5.5\,\text{T}$).
Neutron scattering resolved two sharp excitation modes with mass ratio:

$$\frac{m_2}{m_1} = 1.618 \pm 0.015$$

consistent with the golden ratio $\varphi = 1.618\ldots$

The third particle $m_3$ was also tentatively identified.

## References

- A. B. Zamolodchikov, Int. J. Mod. Phys. A **4**, 4235 (1989) —
  original derivation of E8 integrability.
- G. Mussardo, *Statistical Field Theory*, Oxford University Press
  (2010), Ch. 16 — pedagogical treatment of integrable field theories.
- G. Delfino, J. Phys. A **37**, R45 (2004) — review of Ising field
  theory.
- R. Coldea et al., Science **327**, 177 (2010) — experimental
  observation in CoNb₂O₆.
- V. A. Fateev, Phys. Lett. B **324**, 45 (1994) — E8 Toda S-matrix.

## Used by

- [E8 universality](../universalities/e8.md) — mass ratios stored in QAtlas
- [Ising CFT magnetic perturbation](ising-cft-magnetic-perturbation.md) — physical context
- [TFIM](../models/quantum/tfim.md) — unperturbed theory

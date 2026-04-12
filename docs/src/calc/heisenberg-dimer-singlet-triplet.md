# Heisenberg Dimer: Singlet-Triplet Splitting

## Setup

Two spin-$\tfrac{1}{2}$ particles coupled by the isotropic Heisenberg
exchange interaction:

$$H = J\,\mathbf{S}_1 \cdot \mathbf{S}_2$$

with $J > 0$ (antiferromagnetic) or $J < 0$ (ferromagnetic).
The Hilbert space is $\mathbb{C}^2 \otimes \mathbb{C}^2$ (dimension 4).

## Calculation

### Rewriting in terms of total spin

The total spin is $\mathbf{S}_{\mathrm{tot}} = \mathbf{S}_1 + \mathbf{S}_2$.
Squaring both sides:

$$\mathbf{S}_{\mathrm{tot}}^2
  = \mathbf{S}_1^2 + \mathbf{S}_2^2 + 2\,\mathbf{S}_1 \cdot \mathbf{S}_2$$

so

$$\mathbf{S}_1 \cdot \mathbf{S}_2
  = \frac{1}{2}\bigl[\mathbf{S}_{\mathrm{tot}}^2 - \mathbf{S}_1^2 - \mathbf{S}_2^2\bigr]
  = \frac{1}{2}\bigl[S_{\mathrm{tot}}(S_{\mathrm{tot}}+1) - \tfrac{3}{2}\bigr]$$

where $\mathbf{S}_i^2 = s(s+1) = \tfrac{3}{4}$ for $s = \tfrac{1}{2}$.

### Singlet ($S_{\mathrm{tot}} = 0$)

$$\mathbf{S}_1 \cdot \mathbf{S}_2 = \frac{1}{2}\bigl[0 - \tfrac{3}{2}\bigr]
  = -\frac{3}{4}$$

$$E_s = -\frac{3J}{4}, \qquad
  |\psi_s\rangle = \frac{1}{\sqrt{2}}\bigl(|\!\uparrow\downarrow\rangle
    - |\!\downarrow\uparrow\rangle\bigr)$$

### Triplet ($S_{\mathrm{tot}} = 1$)

$$\mathbf{S}_1 \cdot \mathbf{S}_2 = \frac{1}{2}\bigl[2 - \tfrac{3}{2}\bigr]
  = \frac{1}{4}$$

$$E_t = \frac{J}{4}, \qquad
  |t_m\rangle \in \{|\!\uparrow\uparrow\rangle,\;
  \tfrac{1}{\sqrt{2}}(|\!\uparrow\downarrow\rangle + |\!\downarrow\uparrow\rangle),\;
  |\!\downarrow\downarrow\rangle\}$$

### Gap

$$\Delta = E_t - E_s = \frac{J}{4} - \Bigl(-\frac{3J}{4}\Bigr) = J$$

## Result

$$\boxed{\mathrm{Spectrum} = \left\{-\frac{3J}{4},\;\frac{J}{4},\;\frac{J}{4},\;\frac{J}{4}\right\}}$$

The singlet-triplet gap is $\Delta = J$. For antiferromagnetic
coupling ($J > 0$) the ground state is the singlet; for ferromagnetic
coupling ($J < 0$) the ground state is the triplet.

## References

- Any quantum mechanics textbook, e.g. Sakurai, *Modern Quantum Mechanics*, Ch. 4.
- A. Auerbach, *Interacting Electrons and Quantum Magnetism* (Springer, 1994), Section 2.

## Used by

- [Heisenberg Model](../models/quantum/heisenberg.md)

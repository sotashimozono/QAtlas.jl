# XXZ Chain: Luttinger Parameters from the Bethe Ansatz

## Setup

The spin-$\tfrac{1}{2}$ XXZ chain,

$$H = J \sum_{i} \bigl[\, S^x_i S^x_{i+1} + S^y_i S^y_{i+1}
                       + \Delta\, S^z_i S^z_{i+1} \,\bigr],$$

is Bethe-ansatz integrable for every $\Delta \in \mathbb{R}$.  In the
critical regime $-1 < \Delta \le 1$ we parameterise the anisotropy as

$$\Delta = \cos\gamma, \qquad \gamma \in [0, \pi),$$

so that the isotropic antiferromagnetic point $\Delta = 1$ corresponds
to $\gamma = 0$ and the XX / free-fermion point $\Delta = 0$ to
$\gamma = \pi/2$.

We derive the two Luttinger-liquid parameters

- the compactification parameter $K$, and
- the sound velocity $u$,

both of which are closed-form functions of $\gamma$ alone (up to the
overall $J$ scale).

## Calculation

### Bethe equations and dressed charge

In the $M = N/2$ sector (zero magnetisation) the Bethe rapidities
$\{\lambda_j\}$ satisfy, for $|\Delta| < 1$,

$$\left[\frac{\sinh(\gamma(\lambda_j + i/2))}
              {\sinh(\gamma(\lambda_j - i/2))}\right]^{N}
    = -\prod_{k=1}^{M}
    \frac{\sinh(\gamma(\lambda_j - \lambda_k + i))}
         {\sinh(\gamma(\lambda_j - \lambda_k - i))}.$$

In the thermodynamic limit the ground-state rapidity density
$\rho(\lambda)$ on the real line obeys the linear integral equation

$$\rho(\lambda)
  + \int_{-\infty}^{\infty} K_2(\lambda - \mu)\,\rho(\mu)\,d\mu
  = K_1(\lambda),$$

with kernels

$$K_n(\lambda) = \frac{1}{\pi}\,
    \frac{\sin(n\gamma)}{\cosh(2\gamma\lambda) - \cos(n\gamma)}.$$

Fourier transform under the convention
$\widehat f(k) = \int e^{-ik\lambda} f(\lambda)\,d\lambda$ gives

$$\widehat{K_n}(k) = \frac{\sinh[(\pi - n\gamma)k/(2\gamma)]}
                          {\sinh[\pi k/(2\gamma)]},$$

and the so-called **dressed charge** function $Z(\lambda)$ (the central
object of Luttinger-parameter extraction) obeys

$$Z(\lambda) + \int_{-\infty}^{\infty} K_2(\lambda - \mu)\,Z(\mu)\,d\mu = 1,$$

with $\widehat Z(0) = 1 / [1 + \widehat{K_2}(0)]$.  Using
$\widehat{K_2}(0) = (\pi - 2\gamma)/\gamma$ — which is one of the
classic evaluations in Takahashi (1999), Ch. 4 — gives

$$Z = Z(\infty) = \sqrt{\frac{\pi}{2(\pi - \gamma)}}.$$

### Luttinger parameter

The Luttinger parameter $K$ of a 1D critical chain is the square of
the dressed charge at the Fermi point:

$$\boxed{\;K = Z^{2} = \frac{\pi}{2(\pi - \gamma)}\;}$$

**Check at the canonical points**:

| Point | $\gamma$ | $K$ |
|-------|----------|-----|
| AF Heisenberg, $\Delta = 1$ | $0$ | $\tfrac{1}{2}$ |
| XX, $\Delta = 0$ | $\pi/2$ | $1$ |
| FM boundary, $\Delta \to -1^+$ | $\pi^-$ | $\to \infty$ |

The $\Delta = 0$ value $K = 1$ is the free-fermion (spin-isotropic
bosonisation) result; the $\Delta = 1$ value $K = 1/2$ is the
celebrated "$K = 1/2$ at the SU(2)-symmetric point" needed to explain
the logarithmic correction to the spin-spin correlator at the
Heisenberg point (Affleck, Giamarchi).

### Sound velocity

The Luttinger / spin-wave velocity $u$ is read off from the dressed
energy $\varepsilon(\lambda)$:

$$u = \lim_{\lambda \to \infty}
        \frac{\varepsilon'(\lambda)}{2\pi \rho(\lambda)}.$$

Evaluating with the same Fourier machinery (Takahashi Ch. 4, eq. 4.2.35
or Giamarchi Appendix H) gives the closed form

$$\boxed{\;u(\Delta) = J\cdot \frac{\pi}{2}\,\frac{\sin\gamma}{\gamma}\;}$$

with, as above, $\gamma = \arccos\Delta$.

**Check at the canonical points**:

| Point | $\gamma$ | $u/J$ | Interpretation |
|-------|----------|-------|----------------|
| XX, $\Delta = 0$ | $\pi/2$ | $1$ | Fermi velocity $v_F$ |
| AF Heisenberg, $\Delta = 1$ | $0$ | $\pi/2$ | des Cloizeaux–Pearson |

The $\Delta = 1$ value recovers the classical 1962 des Cloizeaux–Pearson
spinon dispersion $\epsilon(k) = (\pi J/2) |\sin k|$, whose slope at
the Fermi point $k = \pi/2$ is exactly $\pi J/2$.  The $\Delta = 0$
value is the XX-chain tight-binding Fermi velocity $v_F = 2t\sin k_F$
with $t = J/2$ and half-filling $k_F = \pi/2$.

### Smooth limit at $\Delta = 1$

The ratio $\sin\gamma / \gamma$ has a removable singularity at
$\gamma = 0$:

$$\lim_{\gamma \to 0^+} \frac{\sin\gamma}{\gamma} = 1,$$

so $u(\Delta \to 1^-) = (\pi/2) J$ smoothly.  The QAtlas implementation
special-cases $\gamma \approx 0$ to avoid the $0/0$ evaluation in
floating point.

## Result

$$\boxed{\,
K(\Delta) = \frac{\pi}{2(\pi - \gamma)},
\qquad
u(\Delta) = J\cdot \frac{\pi}{2}\,\frac{\sin\gamma}{\gamma},
\qquad
\gamma \equiv \arccos\Delta,
\;-1 < \Delta \le 1.\,}$$

Both are exact, closed-form Luttinger-liquid parameters across the
full critical regime.  Every long-distance exponent of the XXZ chain
(spin-spin correlation decay, dimer correlation decay, finite-size
gap coefficients, …) is a rational function of $K$ alone.

## References

- L. Hulthén, Ark. Mat. Astron. Fys. **26A**, No. 11 (1938) — $\Delta = 1$
  energy.
- J. des Cloizeaux, J. J. Pearson, Phys. Rev. **128**, 2131 (1962) —
  spin-wave dispersion at $\Delta = 1$.
- F. D. M. Haldane, Phys. Rev. Lett. **45**, 1358 (1980);
  Phys. Rev. Lett. **47**, 1840 (1981) — bosonisation of XXZ.
- M. Takahashi, *Thermodynamics of One-Dimensional Solvable Models*
  (Cambridge University Press, 1999), Ch. 4.
- T. Giamarchi, *Quantum Physics in One Dimension* (Oxford, 2004),
  Ch. 6 and Appendix H.
- I. Affleck, Nucl. Phys. B **336**, 517 (1990) — logarithmic
  corrections at the SU(2)-symmetric point.

## Used by

- [XXZ Chain](../models/quantum/xxz.md)
- [Heisenberg Chain](../models/quantum/heisenberg.md) — $\Delta = 1$
  limit reproduces the Hulthén result and the des Cloizeaux–Pearson
  velocity.

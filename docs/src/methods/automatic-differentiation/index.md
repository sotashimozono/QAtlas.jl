# Automatic Differentiation

## What It Is

QAtlas uses forward-mode automatic differentiation (AD) via
ForwardDiff.jl to extract thermodynamic quantities from the
partition function $Z(\beta, J, h, \ldots)$. All equilibrium
observables are derivatives of $\ln Z$:

$$\langle E \rangle = -\frac{\partial \ln Z}{\partial \beta}, \qquad
C_v = \beta^2 \frac{\partial^2 \ln Z}{\partial \beta^2}, \qquad
\langle \sigma\sigma \rangle_{\mathrm{nn}} = \frac{1}{\beta}\frac{\partial \ln Z}{\partial J}$$

By computing $Z$ as a differentiable function and passing dual
numbers through it, all derivatives are obtained exactly (to
machine precision) without finite differences or symbolic
manipulation.

See [calc/ad-thermodynamics-from-z.md](../../calc/ad-thermodynamics-from-z.md)
for the full derivation of thermodynamic identities.

---

## How It Works

ForwardDiff.jl implements forward-mode AD via **dual numbers**.
A dual number $a + b\,\varepsilon$ with $\varepsilon^2 = 0$ carries
a function value and its derivative simultaneously through any
computation:

$$f(a + b\varepsilon) = f(a) + f'(a) \cdot b\varepsilon$$

For second derivatives (specific heat), nested dual numbers are used:
the outer dual differentiates with respect to $\beta$, and the inner
dual differentiates the result again.

```julia
using ForwardDiff

# Energy: first derivative of -ln Z
E(beta) = ForwardDiff.derivative(b -> -log(Z(b)), beta)

# Specific heat: second derivative
Cv(beta) = beta^2 * ForwardDiff.derivative(
    b -> ForwardDiff.derivative(b2 -> log(Z(b2)), b), beta)
```

---

## Key Requirement: `Real` Type Signatures

For AD to work, every function in the computation chain must accept
**any `Real` subtype**, not just `Float64`. The dual number type
`ForwardDiff.Dual{T,V,N} <: Real` must propagate through:

1. Transfer matrix construction
2. Matrix multiplication ($T^{L_x}$)
3. Trace computation
4. Logarithm

Any function with a restrictive type signature
`f(x::Float64)` will fail when called with a `Dual` input.
QAtlas's `fetch` methods are designed to accept `Real` arguments
throughout.

---

## Why `tr(T^Lx)`, Not Eigenvalue Decomposition

The partition function can in principle be computed as
$Z = \sum_i \lambda_i^{L_x}$ from the eigenvalues of $T$. However,
LAPACK's eigenvalue routines (`DSYEV`, `DGEEV`) are implemented in
Fortran and perform internal comparisons and branches that are
incompatible with dual-number arithmetic. Calling `eigvals` on a
matrix of `Dual` elements will fail.

Matrix multiplication and `tr` are purely algebraic operations that
propagate duals correctly. Therefore QAtlas computes

$$Z = \mathrm{tr}(T^{L_x})$$

via repeated squaring ($O(\log L_x)$ matrix multiplications). This
is the reason QAtlas uses the [transfer matrix](../transfer-matrix/index.md)
approach with `tr` rather than eigenvalue decomposition.

---

## Verification

AD-computed derivatives are verified against brute-force ensemble
averages (explicit sum over all $2^N$ configurations) to machine
precision ($\sim 10^{-15}$ relative error). This confirms that:

1. The transfer-matrix construction preserves differentiability.
2. The AD chain rule propagates correctly through `tr`, matrix
   power, and `log`.
3. No non-differentiable operations (e.g., `if` branches on values)
   break the derivative.

---

## Applications in QAtlas

| Observable | Formula | Used in |
|------------|---------|---------|
| Internal energy $\langle E \rangle$ | $-\partial_\beta \ln Z$ | IsingSquare thermal verification |
| Specific heat $C_v$ | $\beta^2 \partial_\beta^2 \ln Z$ | Cross-check #7: $\alpha = 0$ |
| NN correlation $\langle\sigma\sigma\rangle$ | $\beta^{-1}\partial_J \ln Z$ | IsingSquare correlation check |

---

## References

- J. Revels, M. Lubin, T. Papamarkou, "Forward-mode automatic
  differentiation in Julia", arXiv:1607.07892 (2016) ---
  ForwardDiff.jl.
- A. G. Baydin, B. A. Pearlmutter, A. A. Radul, J. M. Siskind,
  "Automatic differentiation in machine learning: a survey",
  JMLR **18**, 1 (2018) --- AD review.

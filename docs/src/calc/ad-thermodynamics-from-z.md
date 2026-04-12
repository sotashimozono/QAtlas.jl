# Thermodynamic Quantities from d(ln Z)/d beta via AD

## Setup

A classical statistical-mechanical system with partition function
$Z(\beta)$ depending on inverse temperature $\beta = 1/T$ and possibly
coupling constants (e.g. $J$, $h$). The partition function may be
computed as $Z = \operatorname{tr}(T^{L_x})$ where $T$ is the
transfer matrix, or from an explicit sum over configurations.

## Calculation

### Key thermodynamic identities

All equilibrium thermodynamic quantities follow from derivatives of
$\ln Z$:

**Internal energy:**

$$\langle E \rangle = -\frac{\partial \ln Z}{\partial \beta}$$

**Specific heat:**

$$C_v = -\beta^2 \frac{\partial^2 \ln Z}{\partial \beta^2}
  = \beta^2\bigl(\langle E^2 \rangle - \langle E \rangle^2\bigr)$$

Equivalently, $C_v = \beta^2 \partial^2(\beta F)/\partial\beta^2$
where $F = -\ln Z / \beta$ is the free energy.

**Nearest-neighbour correlation (Ising):**

$$\langle \sigma_i \sigma_j \rangle_{\mathrm{nn}}
  = \frac{1}{\beta}\frac{\partial \ln Z}{\partial J}$$

**Magnetization in external field $h$:**

$$\langle m \rangle = \frac{1}{\beta}\frac{\partial \ln Z}{\partial h}$$

### Implementation via ForwardDiff.jl

Julia's ForwardDiff.jl implements forward-mode automatic
differentiation (AD) via dual numbers. A dual number
$a + b\,\epsilon$ with $\epsilon^2 = 0$ carries the function value
and its derivative simultaneously through any computation.

For thermodynamic quantities:

```julia
using ForwardDiff

function free_energy(beta, J)
    T = transfer_matrix(beta, J)  # build transfer matrix
    Z = tr(T^Lx)                  # partition function
    return -log(Z) / beta
end

# Energy via first derivative
E(beta) = ForwardDiff.derivative(b -> -log(Z(b)), beta)

# Specific heat via second derivative
Cv(beta) = beta^2 * ForwardDiff.derivative(
    b -> ForwardDiff.derivative(b2 -> log(Z(b2)), b), beta)
```

### Why $\operatorname{tr}(T^{L_x})$, not eigenvalue decomposition

LAPACK's eigenvalue routines (`DGEEV`, `DSYEV`, etc.) do **not**
accept dual-number inputs. The internal Fortran code performs
comparisons and branches that are incompatible with the dual-number
arithmetic. In contrast, matrix multiplication and the trace
operation are purely algebraic and propagate dual numbers correctly.

Therefore, the partition function must be computed as
$Z = \operatorname{tr}(T^{L_x})$ via repeated matrix multiplication,
not as $Z = \sum_i \lambda_i^{L_x}$ from eigenvalues.

For large $L_x$, the matrix power can be computed in $O(\log L_x)$
matrix multiplications via repeated squaring.

### Verification

AD-computed derivatives agree with brute-force ensemble averages
(explicit sum over all $2^N$ configurations) to machine precision
($\sim 10^{-15}$ relative error). This serves as a consistency
check that:

1. The transfer-matrix construction is correct.
2. The AD chain rule is applied correctly through `tr`, matrix
   power, and `log`.
3. No non-differentiable operations (e.g. `if` branches on the
   value) break the derivative propagation.

## Result

$$\boxed{\langle E \rangle = -\frac{\partial \ln Z}{\partial \beta}, \qquad
  C_v = \beta^2\frac{\partial^2 \ln Z}{\partial \beta^2}, \qquad
  \langle \sigma\sigma \rangle_{\mathrm{nn}} = \frac{1}{\beta}\frac{\partial \ln Z}{\partial J}}$$

All thermodynamic observables are obtained from $\ln Z$ by automatic
differentiation. The implementation uses $\operatorname{tr}(T^{L_x})$
rather than eigenvalue decomposition to maintain compatibility with
dual-number arithmetic.

## References

- Standard statistical mechanics; see e.g. K. Huang, *Statistical Mechanics* (Wiley, 1987), Ch. 14.
- J. Revels, M. Lubin, T. Papamarkou, arXiv:1607.07892 (2016) — ForwardDiff.jl.

## Used by

- [Automatic Differentiation Method](../methods/automatic-differentiation/index.md)
- [Ising Square Model](../models/classical/ising-square.md)

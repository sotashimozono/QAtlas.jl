# QAtlas.jl

[![docs: stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://codes.sota-shimozono.com/QAtlas.jl/stable/)
[![docs: dev](https://img.shields.io/badge/docs-dev-purple.svg)](https://codes.sota-shimozono.com/QAtlas.jl/dev/)
[![Julia](https://img.shields.io/badge/julia-v1.12+-9558b2.svg)](https://julialang.org)
[![Code Style: Blue](https://img.shields.io/badge/Code%20Style-Blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

[![codecov](https://codecov.io/gh/sotashimozono/QAtlas.jl/graph/badge.svg?token=dnvvwd7DVo)](https://codecov.io/gh/sotashimozono/QAtlas.jl)
[![Build Status](https://github.com/sotashimozono/QAtlas.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/sotashimozono/QAtlas.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**QAtlas** (QUAntum Reference Table for Exact Tests) is a Julia package providing a curated dictionary of **rigorous results** in quantum and statistical physics. Every stored value is traced to a specific publication and cross-validated against independent calculations.

## Quick Start

```julia
using QAtlas

# Onsager's critical temperature for the 2D Ising model
Tc = QAtlas.fetch(IsingSquare(), CriticalTemperature())  # 2J / ln(1 + √2)

# Yang's spontaneous magnetization
M = QAtlas.fetch(IsingSquare(), SpontaneousMagnetization(); β=0.5)

# 2D Ising universality class: exact critical exponents
e = QAtlas.fetch(Universality(:Ising), CriticalExponents(); d=2)
# (β = 1//8, ν = 1//1, γ = 7//4, η = 1//4, δ = 15//1, α = 0//1, c = 1//2)

# 3D Ising: numerical estimates with uncertainty
e3 = QAtlas.fetch(Universality(:Ising), CriticalExponents(); d=3)
# (α = 0.11009, α_err = 1.0e-5, β = 0.32642, β_err = 1.0e-5, ...)

# Graphene tight-binding spectrum (honeycomb lattice)
λ = QAtlas.fetch(Graphene(), TightBindingSpectrum(); Lx=3, Ly=3, t=1.0)

# Bethe ansatz: Heisenberg chain ground-state energy density
e0 = QAtlas.fetch(Heisenberg1D(), GroundStateEnergyDensity())  # 1/4 − ln 2

# Classical Ising partition function via transfer matrix
Z = QAtlas.fetch(IsingSquare(), PartitionFunction(); Lx=4, Ly=4, β=0.44)
```

## What's Inside

### Models

| Model | Quantities | Method |
| ----- | ---------- | ------ |
| **IsingSquare** | Partition function Z(Lx, Ly, β) | Transfer matrix |
| | Critical temperature T_c | Onsager (1944) exact |
| | Spontaneous magnetization M(T) | Yang (1952) exact |
| **Graphene** (honeycomb TB) | Single-particle spectrum | Bloch Hamiltonian |
| **Kagome** TB | Spectrum + flat band at +2t | Bloch 3×3 |
| **Lieb** TB | Spectrum + flat band at E=0 | Bloch 3×3 |
| **Triangular** TB | Spectrum (frustrated) | Bloch scalar |
| **Heisenberg1D** | Dimer spectrum {−3J/4, J/4³} | Exact |
| | 4-site PBC spectrum | Exact |
| | Ground-state energy density e₀ | Bethe ansatz |
| **TFIM** | Energy, free energy, entropy, ... | BdG + quadrature |

### Universality Classes

Access via `Universality{C}` with dimension `d`:

| Class | d | Type | Key exponents |
| ----- | - | ---- | ------------- |
| **Ising** | 2 | Exact (Rational) | β=1/8, ν=1, c=1/2 |
| | 3 | Numerical (+err) | Conformal bootstrap |
| | ≥4 | Mean-field | Landau |
| **Percolation** | 2 | Exact | β=5/36, ν=4/3 |
| **3-state Potts** | 2 | Exact | β=1/9 |
| **4-state Potts** | 2 | Exact | β=1/12 |
| **XY** | 2 | BKT | η=1/4 |
| | 3 | Numerical | Bootstrap |
| **Heisenberg** | 3 | Numerical | Bootstrap |
| **KPZ** | 1+1D | Exact | β=1/3, z=3/2 |
| **Mean-field** | any | Exact | β=1/2, ν=1/2 |

Exact values use `Rational{Int}` — scaling relations hold with zero floating-point error:

```julia
e = QAtlas.fetch(Universality(:Ising), CriticalExponents(); d=2)
e.α + 2*e.β + e.γ == 2   # Rushbrooke: exactly true
e.γ == e.β * (e.δ - 1)   # Widom: exactly true
```

## Verification

QAtlas is **self-verifying**: every result in `src/` is cross-checked against independent calculations in `test/`. The test suite (800+ tests) includes:

- **Standalone tests**: Special values, scaling relations, limiting cases
- **Lattice verification**: Real-space ED and Bloch diagonalization via [Lattice2D.jl](https://github.com/sotashimozono/Lattice2D.jl) on all 8 supported topologies
- **AD verification**: ForwardDiff-based extraction of thermodynamic quantities from partition functions
- **Cross-verification**: Universality exponents extracted from model-specific results (e.g., β=1/8 from Yang magnetization, z=1 from TFIM gap scaling)

## ForwardDiff Compatibility

The partition function and related functions accept `ForwardDiff.Dual` numbers, enabling automatic differentiation of thermodynamic quantities:

```julia
using ForwardDiff

# Internal energy from partition function
E = -ForwardDiff.derivative(β -> log(QAtlas.fetch(
    IsingSquare(), PartitionFunction(); Lx=3, Ly=3, β=β)), 0.5)

# Heat capacity via second derivative
C = β² * ForwardDiff.derivative(
    β -> ForwardDiff.derivative(β -> log(Z(β)), β), β)
```

## Installation

```julia
using Pkg
Pkg.add("QAtlas")
```

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines on adding new results, writing tests, and the citation standards we follow.

## License

MIT License. See [LICENSE](LICENSE) for details.

# Roadmap and Open Issues

## Remaining Issues

### Physics accuracy of documentation

The `docs/src/calc/` derivation notes provide overview-level arguments
for each result. Several notes — especially around E8 and integrable
field theory — describe **what** is true and **why** it matters, but
do not contain the full calculation. Specific gaps:

- `calc/e8-mass-spectrum-derivation.md`: the S-matrix bootstrap
  equations are stated but not solved. The actual derivation requires
  the CDD factor construction and minimal form factor program
  (Zamolodchikov 1989, Fateev 1994).
- `calc/ising-cft-magnetic-perturbation.md`: Zamolodchikov's proof
  that 8 conserved charges survive the σ perturbation is cited but not
  reproduced. This requires counting integrals of motion order by order
  in perturbation theory.
- `calc/yang-magnetization-toeplitz.md`: the Toeplitz determinant
  evaluation via Szegő's theorem is described at a high level but the
  asymptotic analysis is not shown.

These are candidates for deeper `calc/` sub-notes when expertise is
available.

### Peschel method for σ^z σ^z TFIM

The BdG eigenvectors for the TFIM in the σ^z σ^z convention give
**dual-fermion** correlators (Kramers-Wannier dual), not spin-basis
correlators. Computing the spin-basis entanglement entropy via the
Peschel correlation-matrix method requires proper handling of the
duality boundary conditions.

Current workaround: full ED (O(2^N), N ≤ 16 with 128GB RAM).

Resolution path: implement the duality-aware Peschel method, or switch
to the σ^x σ^x convention for the entanglement-specific code path.

### Central charge extraction precision

Current: c = 1/2 within ~10% at N=14 (TFIM), c = 1 within ~20% at
N=12 (Heisenberg, even-l filtering).

The OBC Calabrese-Cardy formula already includes the log-sin finite-size
correction. Remaining errors come from:

1. Lattice discretization (irrelevant operators, O(1/N))
2. Boundary effects not captured by the strip conformal mapping
3. For Heisenberg: SU(2) alternating corrections (-1)^l f(l)

Potential improvements:
- Use PBC chains (eliminate boundary corrections) — requires PBC ED or
  fixing the Peschel method
- Include subleading corrections in the fit model
- Larger N via sparse Lanczos (ground state only)

### Bond counting convention for small PBC lattices

Lattice2D's `bonds(lat)` double-counts bonds when Lx = 2 or Ly = 2
with PBC (each edge is traversed in both directions). The transfer
matrix and brute-force enumeration use the same convention, so results
agree — but the effective coupling is doubled compared to larger
lattices. This is documented but can surprise new contributors.

### ~~Two API styles~~  ✅ Resolved (v0.13)

All models have been migrated to concrete parameter-carrying structs:

- `TFIM(; J, h)`, `E8()`, `XXZ1D(; J, Δ)`, `Heisenberg1D(; J)`,
  `IsingSquare(; J, Lx, Ly)`, `Honeycomb(; t, Lx, Ly)`,
  `Kagome(; t, Lx, Ly)`, `Lieb(; t, Lx, Ly)`, `Triangular(; t, Lx, Ly)`,
  `Universality{C}`.
- Quantities are concrete structs: `Energy`, `FreeEnergy`,
  `MagnetizationX/Y/Z`, `SusceptibilityXX/YY/ZZ`, `VonNeumannEntropy`,
  `ThermalEntropy`, `RenyiEntropy(α)`, `ZZCorrelation{M}` (parametric),
  `ZZStructureFactor`, `CentralCharge`, `LuttingerParameter`,
  `LuttingerVelocity`, `FermiVelocity`, `E8Spectrum`,
  `TightBindingSpectrum`, `ExactSpectrum`, `GroundStateEnergyDensity`,
  `PartitionFunction`, `CriticalTemperature`, `CriticalExponents`,
  `GrowthExponents`, …
- `const SpinWaveVelocity = LuttingerVelocity` — type-level alias for
  the 1D critical spin chain community.
- BCs `OBC(N)`, `PBC(N)` carry their own size; `Infinite()` unchanged.
- Legacy Symbol-dispatch calls (`fetch(:TFIM, :energy, OBC(); N, J, h)`
  etc.) are routed through the deprecation shims in `src/deprecate/`
  with a `maxlog=3` warning.  The directory can be `git rm`-ed
  wholesale at v1.0 cut.
- Aqua.jl static checks are enabled in `test/test_aqua.jl`
  (ambiguities, deps_compat, stale_deps, piracies).

Version bumped to 0.13.0.  See PRs #40 (foundation), #41 (TFIM + E8),
#42 (classical + tight-binding + Graphene→Honeycomb rename), #44
(XXZ1D).

## Future Directions

### Near-term (infrastructure exists)

**Sparse ED for larger N** (N = 18–22):
- `build_tfim` and `build_spinhalf_heisenberg` currently produce dense
  matrices. Sparse construction + iterative Lanczos (Arpack.jl or
  KrylovKit.jl) would enable ground state computation at N ~ 20.
- This would improve central charge extraction significantly
  (c precision scales as 1/N).
- Estimated effort: small (sparse matrix fill + `eigsolve` call).

**More universality classes from Wikipedia table**:
- Directed percolation (d=1,2,3: numerical, d≥4: MF)
- Self-avoiding walk (d=2: exact, d=3: numerical)
- Conserved directed percolation (Manna class)
- Tricritical Ising (c = 7/10)
- These are all data-entry tasks with the `Universality{C}` framework.

**Finite-size scaling util** (`test/util/fss.jl`):
- Systematic extraction of critical exponents from ED data at multiple
  N values. Binder cumulant crossing, data collapse, etc.
- Would make the universality cross-checks more rigorous.

### Medium-term (requires new packages or significant work)

**ITensors.jl DMRG benchmarks**:
- Compare DMRG ground state energy against QAtlas exact values
  (Heisenberg Bethe ansatz, TFIM BdG).
- Requires ITensors.jl as a test dependency.
- Goal: QAtlas as the "ground truth" for tensor network benchmarks.

**Disordered systems — quantitative IRFP**:
- Extract c_eff = ln 2 from the random TFIM at the Fisher IRFP
  with proper disorder averaging (100+ realizations, N ~ 14).
- Random singlet phase of the Heisenberg chain: verify the
  log-averaged correlation function exponent.
- Harris criterion test: compare clean α to disorder relevance.

**Anderson localization** (1D TB + random on-site potential):
- All states localized in 1D (exact result).
- Localization length ξ ∝ 1/W² for weak disorder.
- Can use existing `build_tight_binding` + random diagonal.

**Entanglement entropy for free fermions** (Peschel done right):
- Implement the correlation-matrix method for the σ^x σ^x convention
  TFIM (direct JW, no KW duality needed).
- O(N³) computation enables N ~ 1000 for precise c extraction.

### Long-term (research-level)

**TFIM + longitudinal field (TFIML)**:
- Numerical verification of E8 mass ratios on finite chains.
- Requires computing the excitation spectrum (not just ground state)
  of the TFIML Hamiltonian and identifying the 8 stable particles.
- Connection to `src/universalities/E8.jl`.

**2D models on Lattice2D**:
- Classical Monte Carlo verification of Onsager T_c, Yang M(T) via
  Lattice2DMonteCarlo.
- Finite-size scaling of the Binder cumulant crossing → extract ν.
- Requires fixing the Lattice2DMonteCarlo LatticeCore compat issue.

**Bethe ansatz beyond ground state**:
- XXZ chain Luttinger parameters (K, u) — **done in v0.13** for the
  full critical regime −1 < Δ ≤ 1.  Ground-state energy density
  currently exposes only the three exact points Δ ∈ {−1, 0, 1}.
- General-Δ Yang-Yang ground-state integral (Δ ∈ (−1, 1) smoothly
  between the exact points) — deferred to a follow-up PR; several
  equivalent forms in the literature differ by nontrivial
  substitutions and sign conventions, requires careful boundary-limit
  verification.
- Heisenberg: des Cloizeaux-Pearson dispersion ε(k) = (π/2)|sin k| —
  needed for dynamical structure factor; `LuttingerVelocity` already
  returns the right Δ = 1 limit.
- Finite-temperature Bethe ansatz (TBA) for thermodynamic quantities.

**Kitaev models**:
- Kitaev chain (1D): topological phase boundary, Majorana edge modes.
  Same BdG framework as TFIM.
- Kitaev honeycomb (2D): exact Majorana solution, A/B phase diagram.
  Requires 2D lattice Majorana implementation.

**Detailed E8 derivation**:
- Fill in the missing calculation steps in `calc/e8-mass-spectrum-derivation.md`:
  bootstrap equation solution, CDD factors, minimal form factors.
- Requires expertise in integrable field theory.

## Development Practices

- **Version bump**: patch per PR, minor for significant features.
- **Testing**: 3-tier (standalone, verification, cross-verification).
  Every `src/` value must have independent test verification.
- **Documentation**: Zettelkasten `calc/` notes for derivations, linked
  from model/method/universality pages.
- **Citations**: each stored value traces to a specific paper + equation.
- **BLAS threading**: `runtests.jl` sets `BLAS.set_num_threads(Sys.CPU_THREADS)`.
- **Julia threads**: invoke `Pkg.test` with
  `julia_args=\`-t auto --heap-size-hint=96G\`` on the 36-core / 128GB
  host to avoid swap thrash; default 1-thread runs OOM on the
  entanglement-central-charge tests.
- **Adaptive sizing**: ED tests use larger N when RAM > 50GB.

module QAtlas

export AbstractModel, Model
export BoundaryCondition, Infinite, PBC, OBC
export AbstractQuantity, Quantity
export fetch

# --- Classical Models ---
export IsingSquare, PartitionFunction, CriticalTemperature, SpontaneousMagnetization

# --- Quantum Models ---
export TFIM                                             # v0.13 concrete struct
export E8                                               # v0.13 concrete struct
export Honeycomb, TightBindingSpectrum                  # v0.13 rename: was Graphene
# NOTE: `Kagome`, `Lieb`, `Triangular` are NOT exported — they conflict
# with Lattice2D's topology types of the same name. Access them as
# `QAtlas.Kagome()` / `QAtlas.Lieb()` / `QAtlas.Triangular()` in code
# that also uses `Lattice2D`.  `Graphene` is kept as a top-level
# backward-compat alias for `Honeycomb` (see src/deprecate/).
export Heisenberg1D, ExactSpectrum, GroundStateEnergyDensity

# --- Core Implementation ---
include("core/alias.jl")
include("core/type.jl")
include("core/quantities.jl")
include("core/pfaffian.jl")

# --- Quantity struct exports (new, axis-explicit naming) ---
export Energy, FreeEnergy, SpecificHeat, MassGap, FidelitySusceptibility
export ThermalEntropy, VonNeumannEntropy, RenyiEntropy
export MagnetizationX, MagnetizationY, MagnetizationZ
export MagnetizationXLocal, MagnetizationZLocal, EnergyLocal
export SusceptibilityXX, SusceptibilityYY, SusceptibilityZZ
export XXCorrelation, YYCorrelation, ZZCorrelation
export XXStructureFactor, YYStructureFactor, ZZStructureFactor
export CentralCharge, LuttingerParameter
export FermiVelocity, LuttingerVelocity, SpinWaveVelocity
export E8Spectrum

# --- Universality Classes ---
export Universality, CriticalExponents, GrowthExponents
export Ising2D, KPZ1D, MeanField  # backward-compatible aliases
include("universalities/Universality.jl")
include("universalities/E8.jl")
include("universalities/MeanField.jl")
include("universalities/Ising2D.jl")
include("universalities/KPZ.jl")
include("universalities/Percolation.jl")
include("universalities/Potts.jl")
include("universalities/ONModel.jl")

# --- Models ---
include("models/TFIM.jl")
include("models/TFIM_dynamics.jl")
include("models/TFIM_thermal.jl")
include("models/TFIM_local.jl")
include("models/classical/IsingSquare.jl")
include("models/quantum/tightbinding/Honeycomb.jl")
include("models/quantum/tightbinding/Kagome.jl")
include("models/quantum/tightbinding/Lieb.jl")
include("models/quantum/tightbinding/Triangular.jl")
include("models/quantum/Heisenberg.jl")

# --- Deprecation shims (legacy API) ---
# Loaded last so they can route into any already-registered concrete
# `fetch` method.  See src/deprecate/README.md.
include("deprecate/legacy_fetch.jl")
include("deprecate/legacy_tfim.jl")
include("deprecate/legacy_e8.jl")
include("deprecate/legacy_honeycomb.jl")
export Graphene                                         # backward-compat alias

end # module QAtlas

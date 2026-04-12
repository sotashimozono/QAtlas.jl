module QAtlas

export AbstractModel, Model
export BoundaryCondition, Infinite, PBC, OBC
export AbstractQuantity, Quantity
export fetch

# --- Classical Models ---
export IsingSquare, PartitionFunction, CriticalTemperature, SpontaneousMagnetization

# --- Quantum Models ---
export Graphene, TightBindingSpectrum
# NOTE: `Kagome`, `Lieb`, `Triangular` are NOT exported — they conflict
# with Lattice2D's topology types of the same name. Access them as
# `QAtlas.Kagome()` / `QAtlas.Lieb()` / `QAtlas.Triangular()` in code
# that also uses `Lattice2D`.
export Heisenberg1D, ExactSpectrum

# --- Core Implementation ---
include("core/alias.jl")
include("core/type.jl")
include("core/pfaffian.jl")

# --- Universality Classes ---
export Ising2D, KPZ1D, MeanField, CriticalExponents
include("universalities/E8.jl")
include("universalities/Ising2D.jl")
include("universalities/KPZ.jl")
include("universalities/MeanField.jl")

# --- Models ---
include("models/TFIM.jl")
include("models/TFIM_dynamics.jl")
include("models/TFIM_thermal.jl")
include("models/classical/IsingSquare.jl")
include("models/quantum/tightbinding/Graphene.jl")
include("models/quantum/tightbinding/Kagome.jl")
include("models/quantum/tightbinding/Lieb.jl")
include("models/quantum/tightbinding/Triangular.jl")
include("models/quantum/Heisenberg.jl")

end # module QAtlas

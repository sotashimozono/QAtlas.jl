module QAtlas

export AbstractModel, Model
export BoundaryCondition, Infinite, PBC, OBC
export AbstractQuantity, Quantity
export fetch

# --- Classical Models ---
export IsingSquare, PartitionFunction

# --- Quantum Models ---
export Graphene, TightBindingSpectrum

# --- Core Implementation ---
include("core/alias.jl")
include("core/type.jl")
include("core/pfaffian.jl")

# --- Universality Classes ---
include("universalities/E8.jl")

# --- Models ---
include("models/TFIM.jl")
include("models/TFIM_dynamics.jl")
include("models/TFIM_thermal.jl")
include("models/classical/IsingSquare.jl")
include("models/quantum/Graphene.jl")

end # module QAtlas

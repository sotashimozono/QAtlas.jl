
"""
    AbstractModel
Struct for defining the model.
It should contain all the necessary information to construct the Hamiltonian and perform calculations,
such as the lattice geometry, interaction parameters, and any other relevant details.
"""
abstract type AbstractModel end

struct Model{M} <: AbstractModel
    params::Dict{Symbol,Any}
end
function Model(name::Symbol; kwargs...)
    canon = canonicalize_model(Val(name))
    return Model{canon}(Dict{Symbol,Any}(kwargs))
end
"""
    AbstractBoundaryCondition
Struct for defining the boundary conditions of the system. Available options include:
- `Infinite`: Infinite system with no boundaries.
- `PBC`: periodic boundary conditions
- `OBC`: open boundary conditions
"""
abstract type BoundaryCondition end
struct Infinite <: BoundaryCondition end
struct PBC <: BoundaryCondition end
struct OBC <: BoundaryCondition end

"""
    AbstractQuantity
"""
abstract type AbstractQuantity end

struct Quantity{Q} <: AbstractQuantity end
function Quantity(q::Symbol)
    canon = canonicalize_quantity(Val(q))
    return Quantity{canon}()
end
Quantity(q::AbstractString) = Quantity(Symbol(q))

"""
    function fetch(model::AbstractModel, quantity::AbstractQuantity; kwargs...)
Fetch the value of a quantity for a given model and boundary condition.
each implementation of `fetch` should be specific to the model and quantity being calculated, and may require additional parameters passed through `kwargs`.
"""
function fetch(
    model::AbstractModel, quantity::AbstractQuantity, bc::BoundaryCondition; kwargs...
)
    return error(
        "QAtlas: No implementation found for Model{$(typeof(model).parameters[1])}, Quantity{$(typeof(quantity).parameters[1])} under $(typeof(bc)) condition.",
    )
end
function fetch(m::Symbol, q::Symbol, bc::BoundaryCondition=Infinite(); kwargs...)
    # Forward kwargs to both Model construction and the fetch dispatch:
    # `Model` stores them all in `params`, and the concrete `fetch`
    # method is free to read either from `params` or from kwargs.
    return fetch(Model(m; kwargs...), Quantity(q), bc; kwargs...)
end

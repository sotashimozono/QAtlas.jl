# deprecate/legacy_honeycomb.jl — keep the `Graphene` name alive after
# the v0.13 rename to `Honeycomb` (QAtlas's graphene tight-binding lived
# under the material name, which collides with Lattice2D's topology
# type of the same name; the lattice is named `Honeycomb` going forward).
#
# A plain `const` alias suffices: `typeof(Graphene()) === Honeycomb`, so
# every existing `fetch(Graphene(), ...)` call dispatches through the new
# `Honeycomb` methods without a warning.

"""
    const Graphene = Honeycomb

Backward-compatible alias for [`Honeycomb`](@ref).  New code should
prefer the lattice-named model.
"""
const Graphene = Honeycomb

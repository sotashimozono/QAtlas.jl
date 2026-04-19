# deprecate/legacy_e8.jl — backward-compat shim for the E8 Symbol API.
#
# Routes the pre-v0.13 call:
#   fetch(Model(:E8), Quantity(:E8_spectrum), Infinite())
#   fetch(:E8, :E8_spectrum, Infinite())
# into the concrete-struct method in src/universalities/E8.jl.

# Legacy Model{:E8} is stateless; just forward.
function fetch(::Model{:E8}, ::Quantity{:E8_spectrum}, bc::Infinite; kwargs...)
    return fetch(E8(), E8Spectrum(), bc; kwargs...)
end

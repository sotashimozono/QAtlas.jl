canonicalize_model(::Val{S}) where {S} = S
canonicalize_quantity(::Val{S}) where {S} = S
# src/core/alias.jl

macro register_aliases(target_func, canon, aliases)
    # --- 修正ポイント：中身の Symbol だけを抽出するヘルパー ---
    unwrap(x) = x isa QuoteNode ? x.value : x

    c = unwrap(canon)
    actual_aliases = if aliases isa Expr && aliases.head == :vect
        unwrap.(aliases.args) # 各要素を unwrap
    else
        [unwrap(aliases)]
    end

    output = quote end

    # 1. 正解（c）の登録
    # ここで改めて QuoteNode に包むことで、生成されるコード内では必ずリテラル :symbol になる
    push!(output.args, :($target_func(::Val{$(QuoteNode(c))}) = $(QuoteNode(c))))

    # 2. 各エイリアスの登録
    for al in actual_aliases
        push!(output.args, :($target_func(::Val{$(QuoteNode(al))}) = $(QuoteNode(c))))
    end

    return esc(output)
end

# Quantity Aliases
@register_aliases canonicalize_quantity :energy [:E, :Energy, :e]
@register_aliases canonicalize_quantity :entanglement_entropy [
    :ee, :EE, :S_vN, :EntanglementEntropy
]
@register_aliases canonicalize_quantity :central_charge [:c, :cc, :CentralCharge]
@register_aliases canonicalize_quantity :mass_gap [
    :gap, :Δ, :Delta, :MassGap, :single_particle_gap, :excitation_gap
]
@register_aliases canonicalize_quantity :zz_corr [:ZZ, :zzcorr, :szsz]
@register_aliases canonicalize_quantity :sz_sz_correlation [
    :szsz_correlation, :zz_correlation_rt, :sz_sz_rt
]
@register_aliases canonicalize_quantity :sx_sx_correlation [
    :sxsx_correlation, :xx_correlation_rt, :sx_sx_rt
]
@register_aliases canonicalize_quantity :sz_sz_spreading [
    :szsz_spreading, :zz_spreading, :sz_sz_lightcone
]

# Thermal observables (TFIM_thermal.jl + TFIM_dynamics.jl thermal extensions)
@register_aliases canonicalize_quantity :free_energy [
    :F, :f, :FreeEnergy, :free_energy_density, :helmholtz_free_energy
]
@register_aliases canonicalize_quantity :entropy [:S, :EntropyDensity, :entropy_density]
@register_aliases canonicalize_quantity :specific_heat [
    :Cv, :cv, :SpecificHeat, :heat_capacity
]
@register_aliases canonicalize_quantity :transverse_magnetization [
    :mx, :Mx, :TransverseMagnetization, :sx_expect, :magnetization_x
]
@register_aliases canonicalize_quantity :transverse_susceptibility [
    :chi_xx, :χ_xx, :TransverseSusceptibility, :susceptibility_x
]
@register_aliases canonicalize_quantity :longitudinal_susceptibility [
    :chi_zz, :χ_zz, :LongitudinalSusceptibility, :susceptibility_z, :uniform_susceptibility
]

# v0.13 velocity family — all resolve to the same canonical Symbol
# (`:luttinger_velocity`) at the alias layer.  The legacy shim in
# `src/deprecate/legacy_xxz.jl` then maps the Symbol to
# `LuttingerVelocity()`; `SpinWaveVelocity` is a type-level alias for
# `LuttingerVelocity`.
@register_aliases canonicalize_quantity :luttinger_velocity [
    :v_LL,
    :vLL,
    :u_LL,
    :uLL,
    :sound_velocity,
    :spin_wave_velocity,
    :v_s,
    :vs,
    :SpinWaveVelocity,
    :fermi_velocity,
    :v_F,
    :vF,
    :FermiVelocity,
    :LuttingerVelocity,
]

@register_aliases canonicalize_quantity :luttinger_parameter [
    :K, :K_LL, :LuttingerParameter
]

@register_aliases canonicalize_quantity :ground_state_energy [
    :e0, :E0, :GroundStateEnergyDensity, :ground_state_energy_density
]
@register_aliases canonicalize_quantity :zz_static_thermal [
    :zz_static, :sz_sz_static, :static_zz_thermal
]
@register_aliases canonicalize_quantity :zz_structure_factor [
    :S_zz, :Szz, :static_structure_factor
]

# Site-local (Vector-valued) thermal observables — TFIM_local.jl
@register_aliases canonicalize_quantity :magnetization_x_local [
    :mx_local, :Mx_local, :transverse_magnetization_local, :sx_local
]
@register_aliases canonicalize_quantity :magnetization_z_local [
    :mz_local, :Mz_local, :longitudinal_magnetization_local, :sz_local
]
@register_aliases canonicalize_quantity :energy_local [
    :E_local, :energy_density_local, :local_energy
]

# Model Aliases
@register_aliases canonicalize_model :TFIM [
    :TransverseFieldIsingModel, :transverseFieldIsingModel
]

@register_aliases canonicalize_model :XXZ1D [:XXZ, :xxz, :xxz1d]

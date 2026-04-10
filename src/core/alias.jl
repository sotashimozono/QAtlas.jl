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

# Model Aliases
@register_aliases canonicalize_model :TFIM [
    :TransverseFieldIsingModel, :transverseFieldIsingModel
]

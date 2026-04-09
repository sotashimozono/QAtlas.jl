# API 設計

## 型階層

```
AbstractModel
  └── Model{S}(params::Dict{Symbol,Any})   # S = canonical model name (phantom type)

AbstractQuantity
  └── Quantity{S}()                         # S = canonical quantity name

BoundaryCondition
  ├── Infinite  # 熱力学極限
  ├── PBC       # 周期境界条件（有限サイズ）
  └── OBC       # 開放境界条件（有限サイズ）
```

phantom type によって `fetch` の多重ディスパッチが成立する。
`Model{:TFIM}` と `Model{:KitaevChain}` は異なる型として扱われる。

## fetch インターフェース

3つの呼び出し形式をサポートする。

```julia
# 形式 1: フル指定
fetch(Model(:TFIM; J=1.0, h=0.5, N=24), Quantity(:energy), OBC(); beta=5.0)

# 形式 2: シンボル略記
fetch(:TFIM, :energy, OBC(); J=1.0, h=0.5, N=24, beta=5.0)

# 形式 3: BC 省略（デフォルト = Infinite）
fetch(:TFIM, :energy; J=1.0, h=0.5, beta=5.0)
```

`beta` のような計算パラメータは kwargs で渡す。
`N`, `J`, `h` などのモデルパラメータも kwargs で渡す（形式 2/3 の場合）。

### 戻り値の型

| 物理量の種類 | 戻り値 |
|------------|--------|
| スカラー量（エネルギー、中心荷など） | `Float64` |
| サイト依存量（局所磁化など） | `Vector{Float64}` |
| 質量スペクトルなど | `Vector{Float64}` |
| 複素相関関数 | `Vector{ComplexF64}` |

## 関数群

| グループ | 関数・マクロ | 役割 |
|---------|------------|------|
| 構築 | `Model(sym; kwargs...)` | モデルインスタンス生成 |
| 構築 | `Quantity(sym)` | 物理量インスタンス生成 |
| 取得 | `fetch(model, qty, bc; kwargs...)` | 物理量を取得 |
| エイリアス | `canonicalize_model(::Val{S})` | モデル別名 → 正規名 |
| エイリアス | `canonicalize_quantity(::Val{S})` | 物理量別名 → 正規名 |
| 登録 | `@register_aliases f canon [aliases]` | 新しいエイリアスを登録 |

## 現在登録されているエイリアス

### モデル

| 正規名 | 別名 |
|-------|------|
| `:TFIM` | `:TransverseFieldIsingModel`, `:transverseFieldIsingModel` |

### 物理量

| 正規名 | 別名 |
|-------|------|
| `:energy` | `:E`, `:Energy`, `:e` |
| `:entanglement_entropy` | `:ee`, `:EE`, `:S_vN`, `:EntanglementEntropy` |
| `:central_charge` | `:c`, `:cc`, `:CentralCharge` |
| `:zz_corr` | `:ZZ`, `:zzcorr`, `:szsz` |
| `:E8_spectrum` | `:E8_mass_ratios`, `:E8_masses`, `:mass_ratios`, `:E8_mass_ratio`, `:mass_ratio` |

## 新しいモデルを追加する方法

1. `src/models/NewModel.jl` を作成
2. `src/core/alias.jl` にエイリアスを登録

```julia
@register_aliases canonicalize_model :NewModel [:NewMod, :new_model]
```

3. `fetch` メソッドを実装

```julia
function QAtlas.fetch(
    model::Model{:NewModel},
    qty::Quantity{:energy},
    bc::Infinite;
    beta::Float64,
    kwargs...
)
    J = model.params[:J]
    # ... 計算
    return result
end
```

4. `src/QAtlas.jl` で `include("models/NewModel.jl")` を追加
5. `test/models/test_NewModel.jl` にサニティチェックを追加

## API 設計に関する将来課題

現行の `Model{S}(params::Dict)` 設計は phantom type + Dict で柔軟だが、
パラメータの型検査が実行時にしか行われない。将来的な代替案として
具体的構造体ベースの設計も検討中（[roadmap.md](roadmap.md) 参照）。

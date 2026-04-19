# API 設計（v0.13 版 / dev メモ）

公開ドキュメント側の正本は `docs/src/` にある。ここは実装者向けの
「なぜこの形か」+ 内部ヘルパ関数の地図。

## 型階層

```
AbstractQAtlasModel                          (abstract, 旧 AbstractModel)
  ├── 量子スピン鎖
  │    ├── TFIM(J::Float64, h::Float64)
  │    ├── XXZ1D(J::Float64, Δ::Float64)
  │    └── Heisenberg1D()                    # フィールド無し、J は fetch kwarg
  ├── 量子 Tight-binding（Lattice2D 名前衝突で export なし — QAtlas.Honeycomb 等で参照）
  │    ├── QAtlas.Honeycomb(t::Float64, Lx::Int, Ly::Int)
  │    ├── QAtlas.Kagome(t, Lx, Ly)
  │    ├── QAtlas.Lieb(t, Lx, Ly)
  │    ├── QAtlas.Triangular(t, Lx, Ly)
  │    └── const Graphene = Honeycomb        # v0.13 rename で残した export
  ├── 古典格子
  │    └── IsingSquare(J::Float64, Lx::Int, Ly::Int)
  ├── E8()                                   # 積分可能場の理論（質量スペクトル）
  ├── Universality{C}                        # C = :Ising / :Percolation / :Potts3 / :Potts4
  │                                          #     / :XY / :Heisenberg / :KPZ / :MeanField …
  ├── Ising2D, KPZ1D, MeanField              # Universality を再エクスポートする alias
  │                                          #   const Ising2D  = Universality{:Ising}
  │                                          #   const KPZ1D    = Universality{:KPZ}
  │                                          #   const MeanField = Universality{:MeanField}
  └── Model{S}(params::Dict)                 # v0.12 以前の Symbol 駆動 API
                                             #   src/deprecate/ に隔離、v1.0 で削除

AbstractQuantity
  ├── 熱力学
  │    └── Energy, FreeEnergy, ThermalEntropy, SpecificHeat, MassGap,
  │        FidelitySusceptibility, PartitionFunction, CriticalTemperature
  ├── 量子情報エントロピー
  │    └── VonNeumannEntropy, RenyiEntropy(α)
  ├── 磁化
  │    └── MagnetizationX / MagnetizationY / MagnetizationZ,
  │        MagnetizationXLocal / MagnetizationZLocal, EnergyLocal,
  │        SpontaneousMagnetization               # IsingSquare 専用（Yang M(T)）
  ├── 帯磁率
  │    └── SusceptibilityXX / SusceptibilityYY / SusceptibilityZZ
  ├── 相関関数（実空間 / 運動量空間）
  │    ├── XXCorrelation{M}, YYCorrelation{M}, ZZCorrelation{M}
  │    │     M = :static | :dynamic | :lightcone | :connected
  │    │     ZZCorrelation(; mode = :static) / ZZCorrelation{:dynamic}() 両形あり
  │    └── XXStructureFactor, YYStructureFactor, ZZStructureFactor
  ├── CFT / Luttinger
  │    ├── CentralCharge, LuttingerParameter
  │    ├── FermiVelocity, LuttingerVelocity
  │    └── const SpinWaveVelocity = LuttingerVelocity   (真の型レベル alias)
  ├── スペクトル
  │    └── E8Spectrum, TightBindingSpectrum, ExactSpectrum,
  │        GroundStateEnergyDensity
  ├── Universality
  │    └── CriticalExponents, GrowthExponents
  └── Quantity{S}()                          # 旧 Symbol 駆動 API、deprecate/ のみ

BoundaryCondition
  ├── Infinite                                 # 熱力学極限（サイズ無し）
  ├── OBC(N::Int)                              # 開放境界 + サイズを内包
  └── PBC(N::Int)                              # 周期境界 + サイズを内包
```

`SpontaneousMagnetization` は `MagnetizationZ` の alias ではなく **独立
struct** で、`IsingSquare` の Yang 磁化 `M(T)` を返す専用 dispatch に
使う（Rational 指数 β = 1/8 も経由）。混ぜないこと。

## fetch インターフェース

```julia
# 正準形（concrete struct 3-重ディスパッチ）
fetch(model::AbstractQAtlasModel,
      quantity::AbstractQuantity,
      bc::BoundaryCondition;
      kwargs...)::ReturnType(model, quantity)
```

- モデル identity と物理パラメータは **model struct 側の field** に
  乗せる。thermal kwargs（`beta`, `T`, `ℓ`, `q`, …）だけ kwargs に残す。
- サイズ `N` は `OBC(N)` / `PBC(N)` に内包されるので kwargs から消えた。
  `fetch(...; N=24)` は deprecation shim 経由でのみ受理される。

### 呼び出し例

```julia
# TFIM
fetch(TFIM(; J=1.0, h=1.0), Energy(), OBC(24); beta=Inf)
fetch(TFIM(; J=1.0, h=1.0), MagnetizationX(), OBC(24); beta=5.0)
fetch(TFIM(; J=1.0, h=1.0), ZZCorrelation{:static}(), OBC(24); beta=5.0, i=1, j=10)
fetch(TFIM(; J=1.0, h=1.0), ZZCorrelation(; mode=:dynamic), OBC(24); i=1, j=10, t=0.5)

# XXZ
fetch(XXZ1D(; J=1.0, Δ=0.5), LuttingerParameter(), Infinite())
fetch(XXZ1D(; J=1.0, Δ=1.0), SpinWaveVelocity(), Infinite())

# Tight-binding（Honeycomb / Kagome / Lieb / Triangular は Lattice2D との
# 名前衝突を避けるため unexport — QAtlas. 接頭辞を付けて呼ぶ）
fetch(QAtlas.Honeycomb(; t=1.0, Lx=6, Ly=6), TightBindingSpectrum(), PBC(36))
# Graphene は後方互換 export（`const Graphene = Honeycomb` で型レベル同一）
fetch(Graphene(; t=1.0, Lx=6, Ly=6), TightBindingSpectrum(), PBC(36))

# Classical
fetch(IsingSquare(; J=1.0, Lx=4, Ly=4), PartitionFunction(); β=0.44)
```

### 旧 Symbol API（deprecation shim）

```julia
# 動く（info log、maxlog=1 per (model, quantity) pair）
fetch(:TFIM, :energy, OBC(); N=24, J=1.0, h=1.0, beta=Inf)
fetch(:XXZ, :spin_wave_velocity, Infinite(); J=1.0, Δ=1.0)
```

shim 実装は `src/deprecate/*.jl`:

| ファイル | 役割 |
|----------|------|
| `legacy_fetch.jl` | `fetch(::Symbol, ::Symbol, bc; ...)` エントリ + `@info` 1-shot |
| `legacy_tfim.jl` | `Model{:TFIM}` → `TFIM(; J, h)` 変換 + Symbol quantity map |
| `legacy_e8.jl` | `Model{:E8}` → `E8()` + `Symbol(:E8_masses)` routing |
| `legacy_honeycomb.jl` | `const Graphene = Honeycomb` alias + `Symbol(:Graphene)` |
| `legacy_xxz.jl` | `Model{:XXZ1D}` → `XXZ1D(; J, Δ)` 変換 + velocity 系 alias |

`@info` は `maxlog = 1` で `(model_sym, quantity_sym)` ペアごとに 1 回出る。
テスト側では:

- **incidental callers**（実際の cross-check で legacy 形を使っていた
  もの）は v0.13.5 までに全て concrete-struct API に port 済み
  （PR #73 / #78）。
- **deliberate shim tests**（alias canonicalization などを検証する
  testset）は `@test_logs (:info, r"symbol-dispatch") ...` で包んで
  deprecation info を CI stdout から除去している。

v1.0 cut で `src/deprecate/` ディレクトリごと `git rm -r` すれば卒業
できる（ポートが完了しているので legacy 依存のある実テストは存在
せず、`@test_logs` で包まれた shim testset は PR でまとめて削除する）。

### legacy quantity alias の扱い

`:TransverseMagnetization` / `:TransverseSusceptibility` /
`:LongitudinalSusceptibility` などは **Symbol レベルの alias**
（`src/core/alias.jl` の `canonicalize_quantity`）で、型レベルの
`const ... = ...` ではない。ユーザー向けには canonical な
`MagnetizationX` / `SusceptibilityXX` / `SusceptibilityZZ` を推奨。

`SpontaneousMagnetization` は IsingSquare の Yang 磁化専用の struct
で `MagnetizationZ` とは**別物**（前者は古典 2D Ising の exact
closed-form、後者は量子スピン演算子の期待値）。

## Quantity 命名の設計根拠

### 軸は明示する (MagnetizationX > TransverseMagnetization)

磁化は tensor 量なので、軸不明の形容詞 (`Transverse`, `Longitudinal`)
は曖昧。`MagnetizationX` / `MagnetizationZ` と書けば模型の
Hamiltonian 表式を読むのと同じ感覚で使える。旧名は
Symbol alias でのみサポート（`canonicalize_quantity` 経由）:

- `:TransverseMagnetization` → canonical `:transverse_magnetization`
  → struct `MagnetizationX` （TFIM thermal dispatch）
- `:TransverseSusceptibility` → `:transverse_susceptibility`
  → struct `SusceptibilityXX`
- `:LongitudinalSusceptibility` → `:longitudinal_susceptibility`
  → struct `SusceptibilityZZ`
- `:SpontaneousMagnetization` → canonical `:spontaneous_magnetization`
  → struct `SpontaneousMagnetization`（IsingSquare 専用、
    `MagnetizationZ` と別物）

### Entropy は種類を分ける

熱力学的エントロピー (−∂F/∂T) と量子情報の von Neumann エントロピー
と Rényi エントロピーは 3 つとも「単位が `[bits]`」で表記が揃って
しまうが、計算経路が全く違う。`Entropy` という一語 alias は曖昧を
招くので**禁止**し、以下のいずれかを明示:

- `ThermalEntropy` — `-∂F/∂T`（自由エネルギーの温度微分）
- `VonNeumannEntropy` — reduced density matrix の `-Tr ρ log ρ`
- `RenyiEntropy(α)` — `(1 - α)⁻¹ log Tr ρ^α`; `α = 1` 極限で vN

### Correlator と Structure factor は分離

実空間 2 点相関 `⟨S^z_i S^z_j⟩` と運動量空間の静的構造因子
`S^{zz}(q)` は、一方は実空間対で他方は Fourier 変換した量。
v0.12 以前は `ZZStructureFactor` と `SzSzCorrelation` が混在して
いたのを整理:

- 実空間: `XXCorrelation{M}`, `YYCorrelation{M}`, `ZZCorrelation{M}`
  パラメータ `M::Symbol`:
  - `:static` — `⟨σ^α_i σ^α_j⟩(β)`
  - `:dynamic` — `⟨σ^α_i(t) σ^α_j(0)⟩`
  - `:lightcone` — `:dynamic` の一種、光円錐上のスプレッディング
  - `:connected` — `⟨σσ⟩ - ⟨σ⟩⟨σ⟩`
- Fourier 空間: `XXStructureFactor`, `ZZStructureFactor`（q, ω 指定）

### Velocity 系は 2 つに分離 + alias

物理起源が違う速度を同名にすると CFT/integrable-model の文脈で
識別不能になるため分ける:

- `FermiVelocity` — 自由フェルミオン / tight-binding / BdG
  （Honeycomb / Kagome / Lieb / Triangular / TFIM の Bogoliubov 速度）
- `LuttingerVelocity` — 1D 臨界相互作用系の bosonization 速度
  （XXZ1D / Heisenberg / TFIM-critical の同一値）
- `const SpinWaveVelocity = LuttingerVelocity` — spin chain コミュニティ
  の慣用名（Haldane / Affleck / des Cloizeaux-Pearson）。型レベルで
  同一なのでディスパッチも struct 同一性も成立。

## Aqua.jl 静的チェック

`test/test_aqua.jl` で実行:

```julia
Aqua.test_all(QAtlas;
    ambiguities        = false,   # 下の "ambiguities (soft)" で報告のみ
    deps_compat        = true,
    stale_deps         = true,
    undocumented_names = false,   # TODO: Infinite / OBC / PBC / Ising2D / KPZ1D /
                                  # Model / Quantity の docstring 整備後に true に戻す
    persistent_tasks   = false,
    piracies           = true,
    unbound_args       = true)

# soft: ambiguities は Aqua.detect_ambiguities(QAtlas; recursive=false) で
# 数だけ報告し、テストは pass させる（stdlib + QuadGK 由来を許容）。
```

`Project.toml` の `[compat]` には全ての `[deps]` および `[extras]`
に対する制約を書く（deps_compat が通る条件）。

README には `[![Aqua QA](...)]` バッジを付けているが、現状は
`undocumented_names = false` なので Aqua.jl 公式基準では partial
compliance。undocumented_names を `true` に戻すのは docstring 補強
PR で別扱い。

## テスト実行

`Pkg.test()` だけだと単一 Julia スレッド + デフォルト ヒープで動く
ので 36-core / 128GB の機で遅い。推奨:

```julia
using Pkg
Pkg.activate(".")
Pkg.test(; julia_args=`-t auto --heap-size-hint=96G`)
```

`runtests.jl` が `BLAS.set_num_threads(Sys.CPU_THREADS)` を実行する
ので、ED の dense eigensolve は常にフル並列。`-t auto` は Julia
レベルの並列ループ（`Threads.@threads`）を有効化する分。

# API 設計（v0.13 版 / dev メモ）

公開ドキュメント側の正本は `docs/src/` にある。ここは実装者向けの
「なぜこの形か」+ 内部ヘルパ関数の地図。

## 型階層

```
AbstractQAtlasModel                          (abstract, 旧 AbstractModel)
  ├── TFIM(J::Float64, h::Float64)           # 量子スピン鎖
  ├── XXZ1D(J::Float64, Δ::Float64)          # 量子スピン鎖
  ├── Heisenberg1D(J::Float64)               # 量子スピン鎖
  ├── E8()                                   # 質量スペクトル
  ├── IsingSquare(J::Float64, Lx::Int, Ly::Int)
  ├── Honeycomb(t::Float64, Lx::Int, Ly::Int) # = Graphene (v0.13 rename)
  ├── QAtlas.Kagome(t, Lx, Ly)                # Lattice2D 名前衝突で unexport
  ├── QAtlas.Lieb(t, Lx, Ly)                  # 同上
  ├── QAtlas.Triangular(t, Lx, Ly)            # 同上
  ├── Universality{C}(…)                     # C = :Ising / :Percolation …
  └── Model{S}(params::Dict) [deprecated]    # v0.12 以前の旧 API

AbstractQuantity
  ├── Energy, FreeEnergy, SpecificHeat, MassGap, FidelitySusceptibility
  ├── ThermalEntropy, VonNeumannEntropy, RenyiEntropy{α}
  ├── MagnetizationX/Y/Z, MagnetizationXLocal/ZLocal, EnergyLocal
  ├── SusceptibilityXX/YY/ZZ
  ├── XXCorrelation{M}, YYCorrelation{M}, ZZCorrelation{M}   (M = :static, :dynamic, :lightcone, :connected)
  ├── XXStructureFactor, YYStructureFactor, ZZStructureFactor
  ├── CentralCharge, LuttingerParameter
  ├── FermiVelocity, LuttingerVelocity
  ├── const SpinWaveVelocity = LuttingerVelocity   (型レベル alias)
  ├── PartitionFunction, CriticalTemperature
  ├── CriticalExponents, GrowthExponents
  ├── E8Spectrum, TightBindingSpectrum, ExactSpectrum
  ├── GroundStateEnergyDensity
  └── Quantity{S}() [deprecated]             # 旧 Symbol 駆動 API

BoundaryCondition
  ├── Infinite                                 # 熱力学極限（サイズ無し）
  ├── OBC(N::Int)                              # 開放境界 + サイズを内包
  └── PBC(N::Int)                              # 周期境界 + サイズを内包
```

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
fetch(TFIM(; J=1.0, h=1.0), ZZCorrelation(:static), OBC(24); beta=5.0)

# XXZ
fetch(XXZ1D(; J=1.0, Δ=0.5), LuttingerParameter(), Infinite())
fetch(XXZ1D(; J=1.0, Δ=1.0), SpinWaveVelocity(), Infinite())

# Tight-binding（Honeycomb は unexport）
fetch(QAtlas.Honeycomb(; t=1.0, Lx=6, Ly=6), TightBindingSpectrum(), PBC(36))
# 旧名も使える（後方互換 const Graphene = Honeycomb）
fetch(Graphene(; t=1.0, Lx=6, Ly=6), TightBindingSpectrum(), PBC(36))

# Classical
fetch(IsingSquare(; J=1.0, Lx=4, Ly=4), PartitionFunction(); β=0.44)
```

### 旧 Symbol API（deprecation shim）

```julia
# 動く（warn maxlog=3）
fetch(:TFIM, :energy, OBC(); N=24, J=1.0, h=1.0, beta=Inf)
fetch(:XXZ, :spin_wave_velocity, Infinite(); J=1.0, Δ=1.0)
```

shim 実装は `src/deprecate/*.jl`:

| ファイル | 役割 |
|----------|------|
| `legacy_fetch.jl` | `fetch(::Symbol, ::Symbol, bc; ...)` エントリ |
| `legacy_tfim.jl` | `Model{:TFIM}` → `TFIM(; J, h)` 変換 + Symbol quantity map |
| `legacy_e8.jl` | `Model{:E8}` → `E8()` + `Symbol(:E8_masses)` routing |
| `legacy_honeycomb.jl` | `const Graphene = Honeycomb` alias + `Symbol(:Graphene)` |
| `legacy_xxz.jl` | `Model{:XXZ1D}` → `XXZ1D(; J, Δ)` 変換 + velocity 系 alias |

v1.0 cut で `src/deprecate/` ディレクトリごと `git rm -r` すれば卒業できる。

## Quantity 命名の設計根拠

### 軸は明示する (MagnetizationX > TransverseMagnetization)

磁化は tensor 量なので、軸不明の形容詞 (`Transverse`, `Longitudinal`)
は曖昧。`MagnetizationX` / `MagnetizationZ` と書けば模型の
Hamiltonian 表式を読むのと同じ感覚で使える。旧名は
deprecation alias として残す:

- `TransverseMagnetization = MagnetizationX`
- `TransverseSusceptibility = SusceptibilityXX`
- `LongitudinalSusceptibility = SusceptibilityZZ`
- `SpontaneousMagnetization = MagnetizationZ`

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
    ambiguities        = false,   # QuadGK / LinearAlgebra との dispatch 多発、別途精査
    deps_compat        = true,
    stale_deps         = true,
    undocumented_names = false,   # TODO: 公開 API の docstring 整備後に true に戻す
    persistent_tasks   = false,
    piracies           = true)
```

`Project.toml` の `[compat]` には全ての `[deps]` および `[extras]`
に対する制約を書く（deps_compat が通る条件）。

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

# XXZ モデル（dev メモ）

## Hamiltonian

```
H = J Σᵢ [ Sˣᵢ Sˣᵢ₊₁ + Sʸᵢ Sʸᵢ₊₁ + Δ Sᶻᵢ Sᶻᵢ₊₁ ]
```

`S = σ/2` 規格。`J > 0` が反強磁性（QAtlas の符号規約）。

パラメータ: `J::Float64` (交換結合, 既定 `1.0`), `Δ::Float64`
(異方性, 既定 `0.0` = XX 点)。

## 相図

| 領域 | 相 | 代表点 / エネルギー |
|------|----|---------------------|
| `Δ < -1` | Gapped FM | — |
| `Δ = -1` | 飽和強磁性 | `e₀/J = -1/4` |
| `-1 < Δ < 1` | Luttinger 液体 (c = 1) | — |
| `Δ = 0` | XX / 自由フェルミオン | `e₀/J = -1/π` |
| `Δ = 1` | SU(2) 対称 AF | `e₀/J = 1/4 - ln 2`（Hulthén 1938）|
| `Δ > 1` | Gapped Néel | — |

## 現在の実装 (v0.13)

| Quantity | BC | 実装状況 | 備考 |
|----------|----|----------|------|
| `Energy` | `Infinite` | **3 点のみ**（Δ ∈ {-1, 0, 1}）| 一般 Δ は NaN + warn |
| `GroundStateEnergyDensity` | `Infinite` | `Energy` のエイリアス | |
| `CentralCharge` | `Infinite` | 閉形式 | `|Δ| < 1` で `1.0`, それ以外 NaN + warn |
| `LuttingerParameter` | `Infinite` | 閉形式 | `K = π / (2(π - γ))`, γ = arccos Δ |
| `LuttingerVelocity` | `Infinite` | 閉形式 | `u = J(π/2) sin γ / γ` |
| `SpinWaveVelocity` | `Infinite` | `LuttingerVelocity` の型エイリアス | `const SpinWaveVelocity = LuttingerVelocity` |

## 保留事項（v0.14 以降）

- **一般 Δ の Yang-Yang 積分**: 文献（Takahashi 1999, Giamarchi 2004,
  原論文 Yang-Yang 1966）で Bethe ansatz 積分表現に複数の等価形がある
  が、符号規約と変数変換が微妙に異なり、Δ = 0 と Δ = 1 の両境界で
  数値一致を取るのに骨が折れる。現在は 3 点を exact に提供し、一般
  Δ は `NaN` + warn を返す実装。
- **有限サイズ補正**: Woynarovich-Eckle (1987) 型の 1/L² 補正は未実装。
- **有限温度自由エネルギー**: Takahashi-Suzuki の Thermal Bethe Ansatz
  (TBA) 方程式は未実装。

## Legacy Symbol API

`src/deprecate/legacy_xxz.jl` 経由で旧 API もサポート:

```julia
QAtlas.fetch(:XXZ, :energy, Infinite(); J=1.0, Δ=0.0)
QAtlas.fetch(:XXZ, :luttinger_parameter, Infinite(); Δ=0.0)
QAtlas.fetch(:XXZ, :spin_wave_velocity, Infinite(); Δ=1.0)
```

対応する velocity 系の Symbol alias: `:v_F`, `:v_LL`, `:u`,
`:fermi_velocity`, `:spin_wave_velocity`, `:sound_velocity`.

## 参考文献

- L. Hulthén, Ark. Mat. Astron. Fys. **26A**, No. 11 (1938) — Δ = 1.
- C. N. Yang, C. P. Yang, Phys. Rev. **150**, 321 (1966) — 一般 Δ の
  Bethe 積分方程式。
- M. Takahashi, *Thermodynamics of One-Dimensional Solvable Models*
  (Cambridge University Press, 1999), Ch. 4.
- T. Giamarchi, *Quantum Physics in One Dimension* (Oxford, 2004),
  Ch. 6 + Appendix H（Luttinger 液体の K, u 公式の出典）。
- F. D. M. Haldane, Phys. Rev. Lett. **45**, 1358 (1980);
  Phys. Rev. Lett. **47**, 1840 (1981) — XXZ の bosonization。

## 公開ドキュメント

- API / 結果: [docs/src/models/quantum/xxz.md]
- 導出: [docs/src/calc/xxz-luttinger-parameters.md]

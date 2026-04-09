# ロードマップ

## 数値データ掲載方針

数値結果を QAtlas に掲載する場合は以下を必須とする:

1. 査読済み論文 or well-known preprint（arXiv）の明示的な reference
2. 測定量・値・不確かさ・測定条件のセット記述
3. 実験値か数値計算かの区別を明記

掲載例（将来）:
```
CoNb₂O₆ での E8 質量比 m₂/m₁
  値: 1.619 ± 0.015（実験）
  条件: 横磁場 Bc = 5.5T 近傍
  出典: Coldea et al., Science 327, 177 (2010)
```

## 将来追加するモデル

現在の主な関心モデルとその既知の厳密解:

| モデル | 解法 | 実装可能な物理量 | 優先度 |
|--------|------|----------------|--------|
| **KitaevChain** | JW + BdG（TFIM と同構造） | エネルギー、ギャップ、位相相境界、エッジモード | 高 |
| **Heisenberg (XXX)** | Bethe ansatz | 基底状態エネルギー密度、スピン波速度 | 高 |
| **XXZ** | Bethe ansatz | エネルギー、Luttinger 液体パラメータ | 中 |
| **KitaevHoneycomb** | Kitaev (2006) の解析的解法 | 基底状態相図、ギャップ | 中 |
| **TFIML (TFIM+L)** | QFT 摂動論（E8 との接続） | 質量ギャップの λ 依存性 | 中 |
| **自由フェルミオン鎖** | 解析的 | 単粒子スペクトル、エンタングルメント | 低 |

### KitaevChain の詳細

```
H = -μ Σᵢ (cᵢ†cᵢ - 1/2) - t Σᵢ (cᵢ†cᵢ₊₁ + h.c.) - Δ Σᵢ (cᵢcᵢ₊₁ + h.c.)
```

- TFIM と同じ JW + BdG 解法 → TFIM の実装を流用可能
- トポロジカル相（\|μ\| < 2t, Δ≠0）でマヨラナエッジモードが出現
- 位相相境界: \|μ\| = 2t（Δ≠0 のとき）
- 参考: Kitaev, *Phys. Usp.* 44, 131 (2001)

### KitaevHoneycomb の詳細

```
H = -Jₓ Σ_{x-links} σˣᵢσˣⱼ - Jᵧ Σ_{y-links} σʸᵢσʸⱼ - Jᵤ Σ_{z-links} σᶻᵢσᶻⱼ
```

- Majorana フェルミオン表示で厳密に解ける（2D）
- 相図: A 相（トポロジカル）/ B 相（ギャップレス）
- 参考: Kitaev, *Ann. Phys.* 321, 2 (2006)

## 将来のサブモジュール案

| モジュール | 内容 | 優先度 | 備考 |
|-----------|------|--------|------|
| **ExactDiag.jl** | 小 N の厳密対角化（参照値として） | 中 | QAtlas のテストに使える |
| **QMC.jl** | 量子モンテカルロの数値結果 | 低 | reference 要件が高い |

## API 設計の再検討

現行設計: `Model{S}(params::Dict{Symbol,Any})`

```julia
fetch(:TFIM, :energy; J=1.0, h=0.5, beta=2.0)
```

**検討中の代替案: 具体的構造体ベース**

```julia
# 各モデルが固有の struct を持つ
struct TFIM <: AbstractModel
    J::Float64
    h::Float64
end

fetch(TFIM(J=1.0, h=0.5), :energy; beta=2.0)
```

**メリット:**
- パラメータの型チェックがコンパイル時に行われる
- FiniteTemperature.jl の `TFIM`, `TFIML` との統一が容易
- IDE の補完が効く

**デメリット:**
- モデルごとに struct 定義が必要
- 現行コードとの互換性が失われる

**TOML パラメータファイル連携（長期課題）:**
configs/*.toml から直接 Model オブジェクトを生成できるようにすると、
FiniteTemperature.jl のパイプラインと QAtlas の参照値比較が自動化できる。

```julia
# 将来のイメージ
spec = load_p1_spec("configs/phase1.toml")
e_ref = fetch(spec.model, :energy; beta=spec.beta_max)
```

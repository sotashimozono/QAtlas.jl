# Ising CFT (c = 1/2)

## 概要

2次元古典 Ising 模型の臨界点 / 1+1次元量子 Ising モデルの量子臨界点は
Virasoro 代数の中心電荷 c=1/2 の共形場理論（Ising CFT）で記述される。

TFIM との対応: `h = J`（臨界点）。

## プライマリ演算子

| 演算子 | 共形次元 h | 物理的対応 |
|--------|-----------|----------|
| 1 (恒等演算子) | 0 | — |
| σ (スピン場) | 1/16 | σᶻ の連続極限 |
| ε (エネルギー密度) | 1/2 | σᶻσᶻ の連続極限 |

スケーリング次元 Δ = h + h̄ = 2h（ユニタリ最小モデル M(3,4)）。

## 臨界指数

| 指数 | 値 | 関係式 |
|------|-----|--------|
| η | 1/4 | 2-η = 2-2Δ_σ |
| ν | 1 | ξ ~ \|g-gc\|^(-ν) |
| β | 1/8 | <σ> ~ (gc-g)^β |
| δ | 15 | <σ> ~ h^(1/δ) at g=gc |

## QAtlas での実装状況

```julia
fetch(:TFIM, :central_charge; J=1.0, h=1.0)   # → 0.5
```

現時点では c のみ。将来的に追加予定：

| 物理量 | 正規名 | 実装 |
|--------|--------|------|
| 中心電荷 | `:central_charge` | ✓ |
| スピン場の scaling dimension | `:scaling_dim_sigma` | 予定 |
| エネルギー密度の scaling dimension | `:scaling_dim_epsilon` | 予定 |

## Ising CFT の摂動

臨界点からの摂動により2種類の massive QFT が得られる：

| 摂動 | 有効理論 | 結果 |
|------|---------|------|
| 横磁場（h ≠ J）| free Majorana fermion | mass gap Δ ~ \|h-J\| |
| 縦磁場（λ ≠ 0）| Zamolodchikov QFT | [E8 質量スペクトル](E8.md) |

## 参考文献

- Di Francesco et al., *Conformal Field Theory* (Springer, 1997), Chapter 12
- Ginsparg, *Applied Conformal Field Theory* (Les Houches Lectures, 1988)

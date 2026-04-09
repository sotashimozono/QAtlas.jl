# QAtlas.jl — 設計ドキュメント

量子多体モデルの厳密解・解析的結果を一元管理するカタログパッケージの設計ノート。

## 掲載方針

| カテゴリ | 掲載条件 |
|---------|---------|
| 厳密解・解析的結果 | 積極的に掲載。導出の概要または参考文献を付記する |
| 数値計算結果 | 査読済み論文または well-known preprint の reference とセットでのみ掲載 |
| 未検証の推測 | 掲載しない |

## クイックスタート

```julia
using QAtlas

# TFIM の熱力学極限エネルギー（β=2.0）
fetch(:TFIM, :energy; J=1.0, h=0.5, beta=2.0)

# 有限サイズ OBC
fetch(:TFIM, :energy, OBC(); J=1.0, h=1.0, N=16, beta=3.0)

# E8 質量比
fetch(:E8, :mass_ratios)
```

エイリアスが使えます：

```julia
fetch(:TFIM, :E; ...)        # :energy の別名
fetch(:TFIM, :ee; ...)       # :entanglement_entropy の別名
```

## ページ一覧

### API・設計

- [api.md](api.md) — 型階層・fetch インターフェース・関数群・拡張方法

### モデル別結果カタログ

- [models/TFIM.md](models/TFIM.md) — 横磁場イジングモデル（実装済み）
- [models/XXZ.md](models/XXZ.md) — XXZ モデル（将来）

### 普遍クラス

- [universalities/E8.md](universalities/E8.md) — E8 例外代数の質量スペクトル（実装済み）
- [universalities/IsingCFT.md](universalities/IsingCFT.md) — Ising CFT (c=1/2)

### 開発ガイド

- [testing.md](testing.md) — テスト哲学・example-based testing の方針
- [roadmap.md](roadmap.md) — 将来モデル・数値データ方針・API 再設計案

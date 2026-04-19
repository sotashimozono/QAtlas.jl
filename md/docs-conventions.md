# Docs conventions（dev memo）

このファイルは QAtlas.jl の docs/src/ 執筆における設計原則を定め、
ドキュメント全体のスタイルと深度を一定に保つことを目的とする。
`md/api.md` / `md/roadmap.md` と同じく **日本語・開発メモ** で、
公開ドキュメントからは外れる。

---

## 1. 二系統構成

docs/src/ は以下の **二系統** で構成する:

| 系統 | ディレクトリ | 目的 |
|------|--------------|------|
| 導出 (calc) | `docs/src/calc/` | 数学的・物理的導出。**1 topic 1 file**, step-by-step で完結 |
| 利用 (API / results / narrative) | `docs/src/models/`, `docs/src/universalities/`, `docs/src/methods/`, `docs/src/verification/` | モデル別 / 普遍クラス別 / 手法別の短い要約と API 例 |

**双方向リンク**を必須とする。calc/ note は自身を利用する
model/universality/method ページへ "Used by" の backlink を持ち、
model 側は自身が主張する結果について "See full derivation:
[...](../../calc/...md)" で該当 calc/ note にリンクする。

## 2. calc/ note の section skeleton

```markdown
# <Title>

## Main result
<boxed 最終結果 + どこで成り立つかを 1–2 文で端的に>

## Setup
<Hamiltonian / 変数定義 / 単位系 / 記号規約 / goal>

## Calculation
### Step 1: <名前付きの operation>
### Step 2: ...
### Step n: Limiting-case checks
<臨界点・自由点・古典極限 2–3 点で代入し既知文献値と比較>

## References
<著者, 誌名 vol, ページ (年) — 引用式番号も明記>

## Used by
- [Consumer page title](../path/to/page.md)
```

**Main result** は note 全体の "abstract"。物理学者が式だけ知り
たくて開いたときに、最初の 1 スクロールで boxed 式と「どの
parameter 領域で成り立つか」の 1 文が得られるようにする。
その後に Setup → Calculation → Checks という順で導出本体を展開し、
最後の Checks セクションで Main result と一致することを verify する。

## 3. model / universality / method ページの section skeleton

```markdown
# <Model / Universality / Method name>

## Overview
<Hamiltonian と parameter, phase diagram の要約, 物理的文脈>

## <Quantity 1>
### Statement
<QAtlas が返す数値 / 式 を boxed で記載>

### References
<論文引用>

### QAtlas API
```julia
QAtlas.fetch(Model(; params), Quantity(), BC(...); kwargs...)
# → value
```

### Verification
| Test file | Method | What is checked |

### Full derivation
See [<calc/ note title>](../../calc/<slug>.md).

## <Quantity 2>
...
```

## 4. Depth standard（calc/ note の拘束条件）

QAtlas に登録されている数値は全て、計算過程が calc/ note で
**step-by-step に追える** ことを目標とする。以下のルールを遵守する:

### 4.1 各行は一つの named transformation で前行から導かれる

式を並べるときは、次の一行が前の一行から:

- *algebraic identity* (代数的恒等式)
- *definition* (定義)
- *Fourier transform* (フーリエ変換)
- *contour integral over poles at $\lambda_n = i(n+1/2)$* (具体的な極)
- *spectral decomposition* (スペクトル分解)
- *change of variable $\mu = \tan(\lambda/2)$* (具体的な変数変換)
- *partial fraction* (部分分数)
- *induction on $n$* (帰納法)

のいずれかで得られることを明示する。ラベルを省略してもよいが、
**読者が一意に再構成できる** レベルの情報は必ず残す。

### 4.2 非自明な代数は行単位で展開する

例: Jordan-Wigner で

$$\tau^z_i \tau^z_{i+1} = -(d^\dagger_i - d_i)(d^\dagger_{i+1} + d_{i+1})$$

を結果だけ書くのは禁止。次のように展開する:

1. 両辺の JW 表示を string \$\prod_{j<i} \tau^z_j\$ 込みで書き下す。
2. 隣接 site での string が \$\tau^z_i\$ の 1 因子ぶんだけ異なる
   ことを指摘し、string が \$\tau^z_i\$ に潰れることを示す。
3. \$\tau^z_i = \tau^+_i + \tau^-_i\$ を代入し、± 符号を追う。
4. 4 項の積を展開し、Pauli 代数の制約で 2 項に合流。

ここまで書いて初めて「等号」が正当化される。

### 4.3 積分は評価する

"using Takahashi Ch. 4 Eq. 4.2.35" のような textbook 転記は禁止。
積分を評価する場合は:

- contour を選び(例: 上半面 +半円)、
- 極を全て列挙し (例: \$\lambda_n = i(n + 1/2),\; n \ge 0\$)、
- residue を一つずつ計算し (例: \$\operatorname{Res}_{\lambda_n}
  = -i(-1)^n / \pi\$)、
- 無限和に収束級数公式を適用する (例: \$\sum_{n \ge 0} (-1)^n/(2n+1)^2
  = G\$, Catalan 定数)。

textbook 結果に頼ってよいのは **genuinely 外部の定理**
(Szegő theorem, Cauchy residue theorem, Stone-Weierstrass) のみ。
この場合も一文で定理文を書いてから適用する。

### 4.4 禁止句

以下の phrase を **使用してはいけない**:

- "it can be shown"
- "we omit details"
- "standard calculation gives"
- "as in [textbook]"
- "one can verify"
- "it is easy to see that"
- "it follows immediately"

一行で済むなら一行書く。長いなら長く書く。

### 4.5 Limiting-case checks

Calculation セクションの末尾は必ず **Limiting-case checks**
sub-section とする。Main result の式を 2–3 の canonical parameter
値で評価し、既知の文献値と比較する。例 (XXZ の Luttinger K):

| Point | \$\gamma\$ | \$K\$ | Literature |
|-------|------------|-------|-----------|
| AF Heisenberg, \$\Delta = 1\$ | \$0\$ | \$1/2\$ | Affleck 1990 |
| XX, \$\Delta = 0\$ | \$\pi/2\$ | \$1\$ | free-fermion |
| FM boundary, \$\Delta \to -1^+\$ | \$\pi^-\$ | \$\to \infty\$ | Haldane 1980 |

boxed 式が 3 点全てで文献値に一致して初めて "Main result" は
verified と呼べる。

## 5. 引用様式

### 5.1 論文

```
P. Pfeuty, Ann. Phys. 57, 79 (1970).
```

形式: `著者, 誌名 vol, ページ (年).`

- vol は**太字なし**（inline の plain text で十分）
- ページは始端のみ (70 ではなく 79)
- 複数著者は `A. B. Smith, C. D. Jones`
- 特定式を引用するときは末尾に `eq. (3.6)` を付ける:

  ```
  P. Pfeuty, Ann. Phys. 57, 79 (1970), eq. (3.6).
  ```

### 5.2 書籍

```
M. Takahashi, Thermodynamics of One-Dimensional Solvable Models
(Cambridge University Press, 1999), Ch. 4.
```

形式: `著者, Title (Publisher, Year), Ch./§.`

### 5.3 Preprint

arXiv ID のみ。`arXiv:2401.12345 [cond-mat.stat-mech]`.

### 5.4 年号は必須

**すべての引用に発行年を必ず付ける**。読者が文献検索で辿れる
最小情報として year はマスト。

## 6. 数式記法

| 用途 | 書き方 |
|------|--------|
| display math | `$$...$$` |
| inline math | `$...$` |
| boxed 最終結果 | `\boxed{...}` |
| equation 番号 | `\tag{3.2}` (同一 page 内で cross-ref するときだけ) |
| 章内参照 | 直近の boxed 式を指すなら `by the boxed equation above` で足りる |

MathJax 3 設定は `docs/make.jl` で `packages = ["base", "ams",
"autoload", "physics"]` を有効化してある。`\bra`, `\ket`, `\abs`,
`\norm` などの physics 記号はそのまま使える。

## 7. 言語方針

- `docs/src/` は **英語のみ**。公開ドキュメントであり、海外
  コミュニティからも読めるよう統一する。
- `md/` は **日本語可**。dev memo であり publish されない。

## 8. 双方向リンク

### 8.1 calc/ note → consumer

末尾に必ず:

```markdown
## Used by

- [TFIM model page](../models/quantum/tfim.md)
- [Jordan-Wigner method](../methods/jordan-wigner/index.md)
```

リンク先は **相対パス**。一つ以上必ず置くこと。

### 8.2 consumer → calc/ note

該当 Quantity の Derivation subsection の末尾に:

```markdown
### Full derivation

See [Jordan-Wigner reduction of the TFIM](../../calc/jw-tfim-bdg.md)
for the complete step-by-step derivation.
```

## 9. File-naming

- `calc/<topic-slug>.md`: kebab-case
- topic slug は **数学的対象** であって用途ではない:
  - `bethe-ansatz-heisenberg-e0.md` ← OK (Bethe ansatz の結果)
  - `heisenberg-ground-state-energy.md` ← NG (どの model / どの手法
    かが不明)
- 1 file = 1 topic を遵守。一つの note が二つの独立した命題を
  主張しはじめたら、分割する。

## 10. 新規 calc/ note 追加時のチェックリスト

PR を出す前に以下を手動 verify:

- [ ] `## Main result` が最上部にあり、boxed 式と 1 文 context
- [ ] `## Setup`, `## Calculation`, `## References`, `## Used by` が揃っている
- [ ] Calculation 末尾に Limiting-case checks がある
- [ ] 禁止句 (§4.4) を grep して 0 件:
      `grep -E "it can be shown|we omit|standard calc|as in [A-Z]" docs/src/calc/<file>.md`
- [ ] `docs/make.jl` の "Derivation Notes" に新 note が追加されている
- [ ] 最低 1 つの model/universality/method ページから forward link がある
- [ ] `julia --project=docs -e 'include("docs/make.jl")'` が Error なく通る

## 11. この標準が適用されるファイル

v0.13 時点の既存 16 件の calc/ note は全てこの standard を満たす
よう段階的に rewrite する (md/roadmap.md の Phase 3a-3p)。
Exemplar は `docs/src/calc/jw-tfim-bdg.md` (Phase 2 で rewrite 済み)。
以降の note は exemplar を参照して depth を揃える。

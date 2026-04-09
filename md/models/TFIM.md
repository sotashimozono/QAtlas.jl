# 横磁場イジングモデル (TFIM)

## Hamiltonian

```
H = -J Σᵢ σᶻᵢ σᶻᵢ₊₁ - h Σᵢ σˣᵢ
```

パラメータ: `J` (Ising 結合), `h` (横磁場)。デフォルト `J=1.0`。

## 相図

```
h/J < 1 : 強磁性秩序相  (<σᶻ> ≠ 0)
h/J = 1 : 量子臨界点    (Ising CFT, c=1/2)
h/J > 1 : 量子常磁性相  (<σᶻ> = 0)
```

臨界点での普遍性クラスは [Ising CFT](../universalities/IsingCFT.md) を参照。
TFIM+L（縦磁場あり）での臨界点摂動については [E8](../universalities/E8.md) を参照。

## 解法: Jordan-Wigner + BdG

Jordan-Wigner 変換でスピン → 自由フェルミオン、BdG 行列を対角化して
準粒子エネルギー {Λₙ} を得る。

**有限サイズ OBC:**
BdG 行列（2N×2N）を数値対角化。固有値の正の半分 {Λₙ}ₙ₌₁ᴺ が準粒子エネルギー。

**熱力学極限 (Infinite/PBC):**
運動量空間での分散関係:
```
Λ(k) = 2 sqrt(J² + h² - 2Jh cos k)
```
熱力学量は k 積分（QuadGK による数値積分, rtol=1e-10）。

## 結果カタログ

| 物理量 | BC | N | 実装 | 手法 |
|--------|-----|---|------|------|
| `energy` | `OBC` | 有限 | ✓ | BdG 数値対角化 |
| `energy` | `Infinite` | ∞ | ✓ | k 積分 (QuadGK) |
| `central_charge` | `Infinite` | ∞ | ✓ | 解析値 c=1/2 |
| `specific_heat` | — | — | 予定 | energy の β 微分 |
| `mass_gap` | `OBC` | 有限 | 予定 | BdG の最小 Λₙ |
| `fidelity_susceptibility` | — | — | 予定 | — |

---

### energy（有限サイズ, OBC）

**数式:**

```
<H>(β) = -Σₙ (Λₙ/2) tanh(β Λₙ / 2)
```

β → ∞ 極限（基底状態）では `beta=Inf` を渡せる。

**呼び出し例:**

```julia
# N=16, h=J=1.0（臨界点）、β=3.0
e = fetch(:TFIM, :energy, OBC(); N=16, J=1.0, h=1.0, beta=3.0)

# β のベクトルも渡せる
betas = 0.5:0.5:5.0
es = fetch(:TFIM, :energy, OBC(); N=16, J=1.0, h=1.0, beta=collect(betas))

# 基底状態
e0 = fetch(:TFIM, :energy, OBC(); N=16, J=1.0, h=0.5, beta=Inf)
```

**期待される傾向:**
- β 増大とともに E が基底状態エネルギーへ単調減少
- N→∞ で Infinite 結果に収束

---

### energy（熱力学極限, Infinite）

**数式:**

```
ε(β) = -(1/π) ∫₀^π dk  (Λ(k)/2) tanh(β Λ(k) / 2)
```

**呼び出し例:**

```julia
e_inf = fetch(:TFIM, :energy; J=1.0, h=0.5, beta=2.0)
# Infinite がデフォルト BC なので BC 省略可
```

**期待される数値（参考, J=h=1.0, β=∞）:**
```
ε₀ ≈ -4/π ≈ -1.2732...  （基底状態エネルギー密度）
```

---

### central_charge（臨界点, Infinite）

**解析値:** c = 1/2（Ising CFT）

臨界点 h = J のときのみ有効。h ≠ J では `NaN` を返し警告を出す。

**呼び出し例:**

```julia
c = fetch(:TFIM, :central_charge; J=1.0, h=1.0)   # → 0.5
c = fetch(:TFIM, :central_charge; J=1.0, h=0.8)   # → NaN + warning
```

## 参考文献

- Sachdev, *Quantum Phase Transitions* (Cambridge, 2011), Chapter 4
- Suzuki, *Prog. Theor. Phys.* 46, 1337 (1971) — TFIM の厳密解

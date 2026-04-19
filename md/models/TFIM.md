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
| `magnetization_x_local` | `OBC` | 有限 | ✓ | Majorana covariance (Vector{Float64}, 長さ N) |
| `magnetization_z_local` | `OBC` | 有限 | ✓ | Z₂ により 0 (Vector{Float64}, 長さ N) |
| `energy_local` | `OBC` | 有限 | ✓ | 対称分割 ε_i, Σ ε_i = ⟨H⟩ (Vector{Float64}, 長さ N) |
| `zz_static_thermal` | `OBC` | 有限 | ✓ | Pfaffian ⟨σᶻ_i σᶻ_j⟩ (Matrix{Float64}, N×N) |
| `specific_heat` | `Infinite` / `OBC` | — | ✓ | `Energy` の β 微分 (ForwardDiff) |
| `mass_gap` | `Infinite` | ∞ | ✓ | 解析値 `Δ = 2\|h − J\|` |
| `mass_gap` | `OBC` | 有限 | ✓ | BdG の最小 Λₙ |
| `entanglement_entropy` (`VonNeumannEntropy`) | `OBC` | 有限 | ✓ | Peschel 法 — Majorana 相関行列 Σ_A の対角化 (O(ℓ³)) |
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

### *_local (サイト依存, OBC)

**動機:**
TPQ / random-sampling による finite-T 推定は各サイトの物理量を生成する
(`mag_Xs::Vector{Float64}`, `corrzzs::Vector{Float64}` など)。bulk 平均 (∑_i / N)
で比較すると self-averaging されて境界・バルクの差が消え、サンプル誤差も
潰れてしまう。exact baseline をサイトごとに出しておけば、typicality の
収束を site-resolved に突き合わせられる。

- `magnetization_x_local` — `[⟨σˣ_i⟩_β for i=1:N]`。Majorana 熱共分散の
  `Σ[2i-1, 2i]` から O(N)。
- `magnetization_z_local` — `[⟨σᶻ_i⟩_β for i=1:N]`。Z₂ 対称性で恒等的に 0
  (Gaussian state, 奇数個の Majorana 積の期待値消失)。random sample の
  揺らぎの exact baseline として使える。
- `energy_local` — 局所エネルギー ε_i。ボンドを両端サイトに対称分割:
    ε_i = -(J/2) (⟨σᶻ_{i-1} σᶻ_i⟩ + ⟨σᶻ_i σᶻ_{i+1}⟩) - h ⟨σˣ_i⟩
  (i=1 と i=N で片側のみ)。`sum(energy_local) == ⟨H⟩` を厳密に満たす。

**呼び出し例:**

```julia
mx_i = fetch(:TFIM, :magnetization_x_local, OBC(); N=16, J=1.0, h=1.0, beta=3.0)
mz_i = fetch(:TFIM, :magnetization_z_local, OBC(); N=16, J=1.0, h=1.0, beta=3.0)
ε_i  = fetch(:TFIM, :energy_local,          OBC(); N=16, J=1.0, h=1.0, beta=3.0)

@assert sum(ε_i) ≈ fetch(:TFIM, :energy, OBC(); N=16, J=1.0, h=1.0, beta=3.0)
@assert sum(mx_i)/16 ≈ fetch(:TFIM, :transverse_magnetization, OBC();
                              N=16, J=1.0, h=1.0, beta=3.0)
```

ZZ 相関の N×N 行列は既存の `:zz_static_thermal` で取得可。

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

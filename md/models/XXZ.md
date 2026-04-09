# XXZ モデル（将来実装予定）

## Hamiltonian

```
H = J Σᵢ (σˣᵢ σˣᵢ₊₁ + σʸᵢ σʸᵢ₊₁ + Δ σᶻᵢ σᶻᵢ₊₁)
```

パラメータ: `J` (交換結合), `Δ` (異方性)。

## 相図

```
Δ < -1 : 強磁性秩序相
Δ = -1 : SU(2) 対称点（強磁性 Heisenberg）
-1 < Δ < 1 : Luttinger 液体相（c=1 CFT）
Δ = 1 : SU(2) 対称点（反強磁性 Heisenberg）
Δ > 1 : Néel 秩序相
```

Δ = 0 が XX モデル（自由フェルミオン）。

## 解法: Bethe ansatz

Bethe ansatz 方程式を数値的に解くことで基底状態エネルギー・励起スペクトルが得られる。

**実装予定の物理量:**

| 物理量 | BC | 実装 |
|--------|-----|------|
| `energy` | `Infinite` | 予定 |
| `central_charge` | `Infinite` | 予定（\|Δ\|<1 で c=1） |
| `luttinger_parameter` | `Infinite` | 予定 |

## 参考文献

- Yang & Yang (1966) — Bethe ansatz の原論文
- Takahashi, *Thermodynamics of One-Dimensional Solvable Models* (Cambridge, 1999)

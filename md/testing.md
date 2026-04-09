# テスト方針

## 哲学

QAtlas のテストは**精密なベンチマーク**ではなく、
「物理的に正しそうかを小さいサイズで確認する」example-based テストを基本とする。

目的: 新しい実装を追加したとき、明らかな物理的誤りを即座に検出できること。
対象: 数値的な厳密性より、傾向・符号・既知の解析値との一致。

## テストカテゴリ

### 1. Dispatch テスト

型構築とエイリアス解決が正しく機能するかを確認する。

```julia
@testset "Model construction" begin
    m1 = Model(:TFIM; J=1.0, h=0.5)
    m2 = Model(:TransverseFieldIsingModel; J=1.0, h=0.5)
    @test typeof(m1) == typeof(m2)   # 同じ phantom type
end

@testset "Quantity aliases" begin
    @test Quantity(:E) == Quantity(:energy)
    @test Quantity(:ee) == Quantity(:entanglement_entropy)
end
```

### 2. 物理的 sanity テスト

符号・単調性・境界条件など物理的に自明な性質を確認する。

```julia
@testset "TFIM energy sanity" begin
    # エネルギーは負（強磁性結合）
    e = fetch(:TFIM, :energy, OBC(); N=8, J=1.0, h=0.5, beta=2.0)
    @test e < 0.0

    # β 増大で単調減少（または一定）
    e1 = fetch(:TFIM, :energy, OBC(); N=8, J=1.0, h=0.5, beta=1.0)
    e2 = fetch(:TFIM, :energy, OBC(); N=8, J=1.0, h=0.5, beta=5.0)
    @test e2 <= e1
end
```

### 3. 既知値テスト

解析的に既知の値と一致するかを確認する。

```julia
@testset "E8 golden ratio" begin
    ratios = fetch(:E8, :mass_ratios)
    ϕ = 2 * cos(π / 5)
    @test ratios[2] / ratios[1] ≈ ϕ  rtol=1e-10
end

@testset "Ising central charge" begin
    c = fetch(:TFIM, :central_charge; J=1.0, h=1.0)
    @test c ≈ 0.5
end
```

### 4. 極限値テスト

高温・低温極限での漸近的な振る舞いを確認する。

```julia
@testset "TFIM high-temperature limit" begin
    # β→0 ではエネルギーは 0 に近づく（等確率混合）
    e_hot = fetch(:TFIM, :energy, OBC(); N=8, J=1.0, h=1.0, beta=0.01)
    @test abs(e_hot) < 0.1
end

@testset "TFIM ground state consistency" begin
    # β→∞（OBC 有限）と基底状態エネルギーの一致
    e_gs  = fetch(:TFIM, :energy, OBC(); N=8, J=1.0, h=0.5, beta=Inf)
    e_low = fetch(:TFIM, :energy, OBC(); N=8, J=1.0, h=0.5, beta=50.0)
    @test e_low ≈ e_gs  rtol=1e-3
end
```

## docstring に example を書く方針

実装する `fetch` メソッドには必ず `# Examples` セクションを設ける。

```julia
"""
    fetch(::Model{:TFIM}, ::Quantity{:energy}, ::OBC; N, J, h, beta) -> Float64

BdG を用いて有限サイズ OBC での TFIM 熱エネルギーを計算する。

# Examples
```julia
julia> fetch(:TFIM, :energy, OBC(); N=8, J=1.0, h=1.0, beta=2.0)
-0.9418...
```
"""
```

ドキュメント中の数値例と `test/` の `@testset` は同じ引数を使い、
ドキュメントが常に動く状態であることを保証する。

## テストファイル構成

```
test/
├── runtests.jl
├── core/
│   ├── test_type.jl        # dispatch テスト
│   └── test_alias.jl       # エイリアステスト
├── models/
│   ├── test_TFIM.jl        # TFIM sanity + 既知値
│   └── test_XXZ.jl         # (将来)
└── universalities/
    ├── test_E8.jl           # E8 質量比
    └── test_IsingCFT.jl    # (将来)
```

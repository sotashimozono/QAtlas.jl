using QAtlas
using Documenter
using Downloads

assets_dir = joinpath(@__DIR__, "src", "assets")
mkpath(assets_dir)
favicon_path = joinpath(assets_dir, "favicon.ico")
logo_path = joinpath(assets_dir, "logo.png")

Downloads.download("https://github.com/sotashimozono.png", favicon_path)
Downloads.download("https://github.com/sotashimozono.png", logo_path)

makedocs(;
    sitename="QAtlas.jl",
    format=Documenter.HTML(;
        canonical="https://codes.sota-shimozono.com/QAtlas.jl/stable/",
        prettyurls=get(ENV, "CI", "false") == "true",
        mathengine=MathJax3(
            Dict(
                :tex => Dict(
                    :inlineMath => [["\$", "\$"], ["\\(", "\\)"]],
                    :tags => "ams",
                    :packages => ["base", "ams", "autoload", "physics"],
                ),
            ),
        ),
        assets=["assets/favicon.ico", "assets/custom.css"],
    ),
    modules=[QAtlas],
    pages=[
        "Home" => "index.md",
        "Models" => [
            "models/index.md",
            "Classical" => [
                "models/classical/index.md",
                "Ising Square" => "models/classical/ising-square.md",
            ],
            "Quantum" => [
                "models/quantum/index.md",
                "TFIM" => "models/quantum/tfim.md",
                "Heisenberg" => "models/quantum/heisenberg.md",
                "Tight-Binding" => [
                    "models/quantum/tightbinding/index.md",
                    "Graphene" => "models/quantum/tightbinding/graphene.md",
                    "Kagome" => "models/quantum/tightbinding/kagome.md",
                    "Lieb" => "models/quantum/tightbinding/lieb.md",
                    "Triangular" => "models/quantum/tightbinding/triangular.md",
                ],
            ],
        ],
        "Universality Classes" => [
            "universalities/index.md",
            "Ising" => "universalities/ising.md",
            "Percolation" => "universalities/percolation.md",
            "Potts" => "universalities/potts.md",
            "KPZ" => "universalities/kpz.md",
            "XY / Heisenberg" => "universalities/on-models.md",
            "Mean-Field" => "universalities/mean-field.md",
            "E8" => "universalities/e8.md",
        ],
        "Verification" => [
            "verification/index.md",
            "Cross-Checks" => "verification/cross-checks.md",
            "Entanglement" => "verification/entanglement.md",
            "Disordered" => "verification/disordered.md",
        ],
        "Methods" => [
            "methods/index.md",
            "Physical" => [
                "Transfer Matrix" => "methods/transfer-matrix/index.md",
                "Bloch Hamiltonian" => "methods/bloch-hamiltonian/index.md",
                "Calabrese-Cardy" => "methods/calabrese-cardy/index.md",
            ],
            "Computational" => [
                "Exact Diagonalization" => "methods/exact-diagonalization/index.md",
                "Automatic Differentiation" => "methods/automatic-differentiation/index.md",
            ],
        ],
    ],
)

deploydocs(; repo="github.com/sotashimozono/QAtlas.jl.git", devbranch="main")

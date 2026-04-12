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
                "Jordan-Wigner" => "methods/jordan-wigner/index.md",
                "Calabrese-Cardy" => "methods/calabrese-cardy/index.md",
            ],
            "Computational" => [
                "Exact Diagonalization" => "methods/exact-diagonalization/index.md",
                "Automatic Differentiation" => "methods/automatic-differentiation/index.md",
            ],
        ],
        "Derivation Notes" => [
            "JW → TFIM BdG" => "calc/jw-tfim-bdg.md",
            "Kramers-Wannier Duality" => "calc/kramers-wannier-duality.md",
            "Transfer Matrix Split" => "calc/transfer-matrix-symmetric-split.md",
            "Yang Magnetization" => "calc/yang-magnetization-toeplitz.md",
            "Heisenberg Dimer" => "calc/heisenberg-dimer-singlet-triplet.md",
            "Bethe Ansatz e₀" => "calc/bethe-ansatz-heisenberg-e0.md",
            "Honeycomb Bloch" => "calc/bloch-honeycomb-dispersion.md",
            "Kagome Flat Band" => "calc/bloch-kagome-flat-band.md",
            "Lieb Flat Band" => "calc/bloch-lieb-flat-band.md",
            "Calabrese-Cardy OBC/PBC" => "calc/calabrese-cardy-obc-vs-pbc.md",
            "AD from ln Z" => "calc/ad-thermodynamics-from-z.md",
            "Scaling Relations" => "calc/ising-scaling-relations.md",
        ],
    ],
)

deploydocs(; repo="github.com/sotashimozono/QAtlas.jl.git", devbranch="main")

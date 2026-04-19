ENV["GKSwstype"] = "100"

using QAtlas, Test, LinearAlgebra, Lattice2D, ForwardDiff, Random
using SparseArrays, KrylovKit
using Aqua

# Use all available BLAS threads for dense eigensolves (ED).
# On multi-core machines this dramatically speeds up eigvals/eigen.
const N_BLAS = min(Sys.CPU_THREADS, 64)
BLAS.set_num_threads(N_BLAS)
println("BLAS threads: $(BLAS.get_num_threads()) / $(Sys.CPU_THREADS) cores")

# Default "fast" test profile keeps every ED at N ≤ 12 (fits in a few
# GB of RAM) so PR CI runs in < 2 minutes on a 4-core GitHub runner.
# Set `QATLAS_TEST_FULL=1` to also run N ∈ {14, 16} sweeps (sparse +
# KrylovKit Lanczos, ~30 min on 128 GB / 36-core hardware).  Nightly
# cron only.
const QATLAS_TEST_FULL = get(ENV, "QATLAS_TEST_FULL", "0") != "0"
println("QATLAS_TEST_FULL = $(QATLAS_TEST_FULL)")
const dirs = ["core/", "universalities/", "models/", "standalone/", "verification/"]

const FIG_BASE = joinpath(pkgdir(QAtlas), "docs", "src", "assets")
const PATHS = Dict()
mkpath.(values(PATHS))

# Load all test utilities ONCE to avoid method-overwrite warnings.
include(joinpath(@__DIR__, "util", "classical_partition.jl"))
include(joinpath(@__DIR__, "util", "tight_binding.jl"))
include(joinpath(@__DIR__, "util", "spinhalf_ed.jl"))
include(joinpath(@__DIR__, "util", "sparse_ed.jl"))
include(joinpath(@__DIR__, "util", "bloch.jl"))
include(joinpath(@__DIR__, "util", "tfim_dense_ed.jl"))

@testset "tests" begin
    test_args = copy(ARGS)
    println("Passed arguments ARGS = $(test_args) to tests.")

    # Aqua static QA first — runs in under a second and catches
    # Project.toml drift (stale / missing compat) before the expensive
    # physics tests.
    @testset "test_aqua.jl" begin
        @time include(joinpath(@__DIR__, "test_aqua.jl"))
    end

    @time for dir in dirs
        dirpath = joinpath(@__DIR__, dir)
        println("\nTest $(dirpath)")
        files = sort(
            filter(f -> startswith(f, "test_") && endswith(f, ".jl"), readdir(dirpath))
        )
        if isempty(files)
            println("  No test files found in $(dirpath).")
            @test false
        else
            for f in files
                @testset "$f" begin
                    filepath = joinpath(dirpath, f)
                    @time begin
                        println("  Including $(filepath)")
                        include(filepath)
                    end
                end
            end
        end
    end
end

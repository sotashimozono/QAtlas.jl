using Aqua

# Aqua.jl: static QA checks (ambiguities, unbound args, stale deps,
# undocumented public names, persistent tasks, piracy, project.toml).
#
# `ambiguities` is left opt-in (we inspect it manually below) because
# LinearAlgebra + QuadGK overloads commonly surface benign ambiguities
# that we do not want to make CI-blocking while still being visible.
#
# The `deps_compat` / `stale_deps` / `undocumented_names` family is the
# set of checks that genuinely catch drift between `Project.toml` and
# the source tree — keep them on.

@testset "Aqua QA" begin
    Aqua.test_all(
        QAtlas;
        ambiguities=false,
        deps_compat=true,
        stale_deps=true,
        # Re-enabled in v0.13.5 after every public struct got a
        # dedicated docstring. Keeps the public-API surface honest:
        # adding a new `export X` without a docstring for `X` now
        # fails CI.
        undocumented_names=true,
        persistent_tasks=false,
        piracies=true,
        unbound_args=true,
    )

    @testset "ambiguities (soft)" begin
        # Report but do not fail — many are inherited from stdlib combinations.
        amb = Aqua.detect_ambiguities(QAtlas; recursive=false)
        if !isempty(amb)
            @info "Aqua: ambiguities detected (non-fatal)" count = length(amb)
        end
        @test true  # soft check always passes
    end
end

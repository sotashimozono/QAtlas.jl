using Test

@testset "QAtlas.jl Core Architecture Tests" begin
    # --- 1. test of multiple dispatch ---
    # Dispatch methods accept `kwargs...` so the Symbol-based convenience
    # wrapper `fetch(:Dummy, :Answer; …)` can forward arbitrary keyword
    # arguments without triggering MethodError on extra params.
    function QAtlas.fetch(
        m::QAtlas.Model{:Dummy}, ::QAtlas.Quantity{:Answer}, ::QAtlas.Infinite; kwargs...
    )
        return 42
    end
    function QAtlas.fetch(
        m::QAtlas.Model{:Dummy},
        ::QAtlas.Quantity{:ParamCheck},
        ::QAtlas.Infinite;
        kwargs...,
    )
        return m.params[:val] * 2
    end
    function QAtlas.fetch(
        ::QAtlas.Model{:Dummy}, ::QAtlas.Quantity{:Answer}, ::QAtlas.OBC; kwargs...
    )
        return "Open Boundary"
    end

    # --- 2. verify definitions ---

    @testset "Construction & Parameters" begin
        m = QAtlas.Model(:Dummy; a=1, b="hello")
        @test m isa QAtlas.Model{:Dummy}
        @test m.params[:a] == 1
        @test m.params[:b] == "hello"

        q = QAtlas.Quantity(:Test)
        @test q isa QAtlas.Quantity{:Test}
    end

    @testset "Dispatch Logic (Core API)" begin
        # fetch(:Dummy, :Answer) -> fetch(Model(:Dummy), Quantity(:Answer), Infinite())
        @test QAtlas.fetch(:Dummy, :Answer) == 42
        @test QAtlas.fetch(:Dummy, :Answer, QAtlas.Infinite()) == 42

        # dispatch of bounary conditions
        @test QAtlas.fetch(:Dummy, :Answer, QAtlas.OBC()) == "Open Boundary"
    end

    @testset "Keyword Argument Passing" begin
        # fetch(:Dummy, :ParamCheck; val=10) 
        # -> Model(:Dummy, val=10) が作られ、params[:val] が参照されるか
        @test QAtlas.fetch(:Dummy, :ParamCheck; val=10) == 20
        @test QAtlas.fetch(:Dummy, :ParamCheck; val=100) == 200
    end

    @testset "Error Handling (Strictness)" begin
        # if no method matches, it should throw a MethodError, not a generic error
        @test_throws ErrorException QAtlas.fetch(:Dummy, :NotExist)
        # if model or quantity is unknown, it should also throw a MethodError
        @test_throws ErrorException QAtlas.fetch(:UnknownModel, :Answer)
    end
    #=
    @testset "Type Stability" begin
        @test @inferred(QAtlas.fetch(:Dummy, :Answer)) == 42
    end
    =#
end

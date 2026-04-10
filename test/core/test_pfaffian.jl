using LinearAlgebra: det

@testset "pfaffian" begin
    @testset "edge cases" begin
        @test QAtlas.pfaffian(zeros(0, 0)) == 1.0
        @test QAtlas.pfaffian(zeros(3, 3)) == 0.0      # odd → 0
        @test QAtlas.pfaffian(zeros(4, 4)) == 0.0      # rank-deficient
        @test_throws ArgumentError QAtlas.pfaffian(zeros(2, 3))
    end

    @testset "2×2" begin
        a = 3.7
        @test QAtlas.pfaffian([0.0 a; -a 0.0]) ≈ a
        @test QAtlas.pfaffian([0.0+0im (1.0+2im); -(1.0+2im) 0.0]) ≈ 1.0 + 2im
    end

    @testset "4×4 closed form" begin
        # Pf(A) = A[1,2] A[3,4] - A[1,3] A[2,4] + A[1,4] A[2,3]
        a, b, c, d, e, f = 1.0, 2.0, 3.0, 4.0, 5.0, 6.0
        A = [
            0.0 a b c;
            -a 0.0 d e;
            -b -d 0.0 f;
            -c -e -f 0.0
        ]
        @test QAtlas.pfaffian(A) ≈ a * f - b * e + c * d   # = 6 - 10 + 12 = 8

        # complex variant
        ac = 1.0 + 1im
        bc = 2.0
        cc = 3.0 - 1im
        dc = 4.0 + 0.5im
        ec = 5.0
        fc = 6.0 + 2im
        Ac = ComplexF64[
            0 ac bc cc;
            -ac 0 dc ec;
            -bc -dc 0 fc;
            -cc -ec -fc 0
        ]
        @test QAtlas.pfaffian(Ac) ≈ ac * fc - bc * ec + cc * dc
    end

    @testset "Pf(A)^2 == det(A) for random skew matrices" begin
        # Real
        for n in (2, 4, 6, 8, 10)
            X = randn(n, n)
            A = (X - X') / 2
            pf = QAtlas.pfaffian(A)
            @test pf^2 ≈ det(A) rtol = 1e-10
        end
        # Complex
        for n in (2, 4, 6, 8, 10)
            X = randn(ComplexF64, n, n)
            A = (X - transpose(X)) / 2     # antisymmetric (NOT skew-Hermitian)
            pf = QAtlas.pfaffian(A)
            @test pf^2 ≈ det(A) rtol = 1e-9
        end
    end

    @testset "Pfaffian sign under row+col swap" begin
        # Swapping rows i,j and columns i,j of an antisymmetric matrix
        # multiplies its Pfaffian by -1.
        for n in (4, 6, 8)
            X = randn(n, n)
            A = (X - X') / 2
            B = copy(A)
            i, j = 2, n - 1
            B[[i, j], :] = B[[j, i], :]
            B[:, [i, j]] = B[:, [j, i]]
            @test QAtlas.pfaffian(B) ≈ -QAtlas.pfaffian(A) rtol = 1e-10
        end
    end

    @testset "Pfaffian of block-diagonal" begin
        # Pf(diag(A_1, A_2)) = Pf(A_1) Pf(A_2)
        a = 1.7
        b = -2.3
        A = [0.0 a 0 0; -a 0 0 0; 0 0 0 b; 0 0 -b 0]
        @test QAtlas.pfaffian(A) ≈ a * b
    end
end

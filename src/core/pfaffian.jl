# ---------------------------------------------------------------------
# Pfaffian of an antisymmetric (skew-symmetric) matrix.
#
# Implementation: Parlett–Reid style congruence transform that takes
# A → T A Tᵀ into a block-diagonal skew form (2×2 blocks), from which
#     Pf(A) = sign(permutation) · ∏ₖ A'[2k-1, 2k]
# after the transformation.  The algorithm is O(n³), works for both
# real and complex `T <: Number`, and only requires `A` to satisfy
# `A[i,j] == -A[j,i]` (antisymmetric, *not* skew-Hermitian).
#
# Partial pivoting (largest-|·| in the remaining column) is used for
# numerical stability.  Each row/column swap contributes a factor of
# -1 to the Pfaffian; each leading pivot `A[k, k+1]` is multiplied in
# directly.  See e.g. Wimmer (2012) "Efficient numerical computation
# of the Pfaffian for dense and banded skew-symmetric matrices".
# ---------------------------------------------------------------------

using LinearAlgebra: LinearAlgebra

export pfaffian

"""
    pfaffian(A::AbstractMatrix) -> eltype(A)

Compute the Pfaffian of an antisymmetric matrix `A`.

`A` must be square of even dimension.  Antisymmetry (`A[i,j] ==
-A[j,i]`) is assumed; only the strictly upper (or lower) triangle is
semantically used — diagonal entries are ignored.

For an odd-dimensional antisymmetric matrix the Pfaffian is zero by
definition, and that is returned.

Works for both real and complex numeric element types.  The algorithm
is numerically stable with partial pivoting and runs in O(n³).
"""
function pfaffian(A::AbstractMatrix{T}) where {T<:Number}
    n = size(A, 1)
    size(A, 2) == n || throw(ArgumentError("pfaffian: matrix must be square"))
    isodd(n) && return zero(T)
    n == 0 && return one(T)
    return _pfaffian_parlett_reid!(Matrix{T}(A))
end

function _pfaffian_parlett_reid!(A::Matrix{T}) where {T<:Number}
    n = size(A, 1)
    result = one(T)
    for k in 1:2:(n - 1)
        # Partial pivoting: find the largest |A[i, k]| for i ∈ k+1 … n
        pivot_idx = k + 1
        pmax = abs(A[k + 1, k])
        for i in (k + 2):n
            if abs(A[i, k]) > pmax
                pmax = abs(A[i, k])
                pivot_idx = i
            end
        end

        # Matrix is rank-deficient along this column ⇒ Pfaffian = 0.
        if pmax == 0
            return zero(T)
        end

        # Swap rows and columns k+1 ↔ pivot_idx to bring pivot to (k, k+1)
        if pivot_idx != k + 1
            for j in 1:n
                A[k + 1, j], A[pivot_idx, j] = A[pivot_idx, j], A[k + 1, j]
            end
            for i in 1:n
                A[i, k + 1], A[i, pivot_idx] = A[i, pivot_idx], A[i, k + 1]
            end
            result = -result
        end

        piv = A[k, k + 1]
        result *= piv
        inv_piv = one(T) / piv

        # Eliminate entries A[i, k] for i > k+1 via
        #   row_i ← row_i - factor * row_{k+1}
        #   col_i ← col_i - factor * col_{k+1}    (congruence transform)
        # where factor = A[i, k] / A[k+1, k] = -A[i, k] / A[k, k+1]
        # (using A[k+1, k] = -A[k, k+1]).
        for i in (k + 2):n
            factor = -A[i, k] * inv_piv
            iszero(factor) && continue
            @inbounds for j in 1:n
                A[i, j] -= factor * A[k + 1, j]
            end
            @inbounds for l in 1:n
                A[l, i] -= factor * A[l, k + 1]
            end
        end
    end
    return result
end

# ─────────────────────────────────────────────────────────────────────────────
# test/util/tfim_dense_ed.jl
#
# Shared dense ED helpers for TFIM tests (test_TFIM_dynamics.jl and
# test_TFIM_thermal.jl). Loaded once from runtests.jl to avoid
# method-overwrite warnings.
# ─────────────────────────────────────────────────────────────────────────────

const _SX = ComplexF64[0 1; 1 0]
const _SZ = ComplexF64[1 0; 0 -1]
const _ID = ComplexF64[1 0; 0 1]

function _op_site(o::AbstractMatrix, k::Int, N::Int)
    m = (k == 1) ? o : _ID
    for j in 2:N
        m = kron(m, j == k ? o : _ID)
    end
    return m
end

function _build_tfim_dense(N::Int, J::Float64, h::Float64)
    d = 2^N
    H = zeros(ComplexF64, d, d)
    for i in 1:(N - 1)
        H -= J * _op_site(_SZ, i, N) * _op_site(_SZ, i + 1, N)
    end
    for i in 1:N
        H -= h * _op_site(_SX, i, N)
    end
    return Hermitian(H)
end

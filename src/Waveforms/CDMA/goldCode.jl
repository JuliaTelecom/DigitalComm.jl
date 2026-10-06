# ----------------------------------------------------
# --- Gold code generation (added for #10 discussion)
# ----------------------------------------------------
"""
Table of "preferred pair" LFSR feedback polynomials used to build Gold code
sequences, indexed by register length `n`. Each entry is a tuple of two tap
vectors (polynomial exponents, including the degree-`n` term; the implicit
constant term `+1` is not listed). Preferred pairs only exist for register
lengths where `n` is odd or `n ≡ 2 mod 4` [1].

[1] R. Gold, "Optimal binary sequences for spread spectrum multiplexing,"
    IEEE Trans. Inf. Theory, vol. 13, no. 4, pp. 619-621, Oct. 1967.
"""
const GOLD_PREFERRED_PAIRS = Dict(
    5  => ([5, 2],          [5, 4, 3, 2]),
    6  => ([6, 1],          [6, 5, 2, 1]),
    7  => ([7, 3],          [7, 3, 2, 1]),
    9  => ([9, 4],          [9, 6, 4, 3]),
    10 => ([10, 3],         [10, 9, 8, 6, 3, 2]),  # GPS C/A G1 / G2
    11 => ([11, 2],         [11, 8, 5, 2]),
)

"""Generate a maximal-length (m-)sequence of length `2^n - 1` from a Fibonacci
LFSR with feedback polynomial exponents `taps` (the degree-`n` term included,
the constant `+1` term implicit). Returns a `BitVector`.
"""
function _mseq(n::Int, taps::Vector{Int})
    L = (1 << n) - 1
    state = (1 << n) - 1 # all-ones seed (any nonzero seed gives the same cyclic sequence)
    out = falses(L)
    for i in 1:L
        out[i] = isodd(state)
        fb = 0
        for t in taps
            fb ⊻= (state >> (n - t)) & 1
        end
        state = (state >> 1) | (fb << (n - 1))
    end
    return out
end

""" Gold code generation.
Generates `2^n - 1` Gold codes of length `2^n - 1`, built from a "preferred
pair" of m-sequences `m1`/`m2` as `m1 .⊻ circshift(m2,k)` for `k = 0 .. L-1`.
Unlike OVSF, Gold codes are not perfectly orthogonal (off-peak correlation is
bounded, not zero), but they keep a *bounded* cross-correlation under
arbitrary time shifts, which OVSF does not — useful for asynchronous/multipath
CDMA links. \\
codes = gold(n, taps1=nothing, taps2=nothing) \\
Input
- n : LFSR register length. Code length/count is `2^n - 1` (not required to
  be a power of 2, unlike `ovsf`). A built-in preferred pair is used if `n`
  is in `keys(GOLD_PREFERRED_PAIRS)`; otherwise `taps1`/`taps2` must be given.
- taps1, taps2 : Optional explicit preferred-pair feedback polynomial taps.
Output
- codes : Matrix of size `(2^n-1) x (2^n-1)` with type Int8, entries ∈ {-1,1}
"""
function gold(n::Int, taps1::Union{Nothing,Vector{Int}}=nothing, taps2::Union{Nothing,Vector{Int}}=nothing)
    if isnothing(taps1) || isnothing(taps2)
        haskey(GOLD_PREFERRED_PAIRS, n) || error("No built-in preferred pair for n=$n. Supply taps1/taps2 explicitly, or use one of $(sort(collect(keys(GOLD_PREFERRED_PAIRS))))")
        taps1, taps2 = GOLD_PREFERRED_PAIRS[n]
    end
    L  = (1 << n) - 1
    m1 = _mseq(n, taps1)
    m2 = _mseq(n, taps2)
    C  = zeros(Int8, L, L)
    for k in 0:L-1
        C[:, k+1] .= 1 .- 2 .* Int8.(m1 .⊻ circshift(m2, k))
    end
    return C
end

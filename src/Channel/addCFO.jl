"""
Adding Carrier Frequency Offset δ in Hz to the input signal x sampled at frequency samplingRate (In Hz) with initial phase ϕ (default 0)
This function does not mutate the input signal. See addCFO! for mutating function.
In case you want to add normalized CFO, set the samplingRate to 1.
# --- Syntax 
y = addCFO(x,delta,samplingRate)
# --- Input parameters 
- x : Input signal 
- δ : Carrier frequency offset [Hz]
- samplingRate : Sampling rate [Hz]
# --- Output parameters 
- y : Signal with CFO
""" 
function addCFO(x::Union{AbstractVector{Complex{T}},AbstractVector{T}},δ,samplingRate,ϕ=0) where {T<:Real}
    y = zeros(Complex{T},length(x))
    addCFO!(y,x,δ,samplingRate,ϕ)
    return y
end

function addCFO!(y::AbstractVector{Complex{T}},x::AbstractVector,δ::Number,samplingRate::Number,ϕ=0) where {T<:Real}
    # --- Basic array check 
    @assert length(x) == length(y) "Input and output should have same length (here input x is $(length(x)) and pre-allocated output has size $(length(y))"
    # --- CFO pulsation 
    ω = δ / samplingRate 
    # --- Adding CFO 
    @inbounds @simd for n ∈ eachindex(x)
        y[n] = x[n] * Complex{T}(exp(2im*π*ω*( (n-1) + ϕ)))
    end
end


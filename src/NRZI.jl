"""
	encodeNRZI(bits::AbstractVector, transitions::Symbol=:low)::AbstractVector

Map a bit sequence to Non-Return-to-Zero Inverted (NRZI) encoded bits.
Expects a vector of bits, e.g. `[0, 1, 1, 0, 0, 1, 0, 1, 0, 1]`.

# Arguments

- `bits::AbstractVector`: Vector of bits to encode.
- `transitions::Symbol`: Symbol indicating the symbol to transition on (`:low`/`:high`). Defaults to `:low`.

# Returns

- `encoded_bits::AbstractVector`: Vector of encoded bits.

# Examples

```jldoctest
julia> encoded_bits = encodeNRZI(Int32[0, 1, 1, 0, 0, 1], :low);

julia> transpose(encoded_bits)
1×6 transpose(::Vector{Int32}) with eltype Int32:
 1  1  1  0  1  1
```

The example below shows how the `transitions` argument affects the encoded bit sequence.

```jldoctest
julia> encoded_bits = encodeNRZI(Int32[0, 1, 1, 0, 0, 1], :high);

julia> transpose(encoded_bits)
1×6 transpose(::Vector{Int32}) with eltype Int32:
 0  1  0  0  0  1
```
"""
function encodeNRZI(bits::AbstractVector, transitions::Symbol=:low)::AbstractVector    
    encoded_bits = similar(bits)
    encodeNRZI!(encoded_bits,bits,transitions)
    return encoded_bits
end

function encodeNRZI!(encoded_bits::AbstractVector,bits::AbstractVector,transitions::Symbol=:low)
    @assert size(encoded_bits) == size(bits) "With NRZI encoding, input and output should have same size ($(size(encoded_bits)) ≂̸ $(size(bits))"
    last_bit = 0
    transition_bit = (transitions == :high) ? 1 : 0
    for n ∈ eachindex(bits)
        if bits[n] == transition_bit
            last_bit =  last_bit ⊻ 1
        end
        encoded_bits[n] = last_bit
    end
    return nothing
end

"""
	decodeNRZI(bits::AbstractVector, transitions::Symbol=:low)::AbstractVector

Decode a  Non-Return-to-Zero Inverted (NRZI) encoded bit sequence.
Expects a vector of bits, e.g. `[0, 1, 1, 0, 0, 1, 0, 1, 0, 1]`.

# Arguments

- `bits::AbstractVector`: Vector of bits to encode.
- `transitions::Symbol`: Symbol represented by a transition in the NRZI coded sequence (`:low`/`:high`). Defaults to `:low`.

# Returns

- `decoded_bits::AbstractVector`: Vector of decoded bits. The first bit of the output depends on a value of a memory bit in the decoder.
  this value is set to `0`.

# Examples

```jldoctest
julia> decoded_bits = decodeNRZI(Int32[1, 1, 1, 0, 1, 1], :low);

julia> transpose(decoded_bits)
1×6 transpose(::Vector{Int32}) with eltype Int32:
 0  1  1  0  0  1
```

The example below shows how the `transitions` argument affects the decoded bit sequence.

```jldoctest
julia> decoded_bits = decodeNRZI(Int32[0, 1, 0, 0, 0, 1], :high);

julia> transpose(decoded_bits)
1×6 transpose(::Vector{Int32}) with eltype Int32:
 0  1  1  0  0  1
```
"""
function decodeNRZI(encoded_bits::AbstractVector, transitions::Symbol=:low)::AbstractVector
    decoded_bits = similar(encoded_bits)
    decodeNRZI!(decoded_bits,encoded_bits,transitions)
    return decoded_bits
end


function decodeNRZI!(decoded_bits::AbstractVector,encoded_bits::AbstractVector, transitions::Symbol=:low)
    @assert size(encoded_bits) == size(decoded_bits) "With NRZI encoding, input and output should have same size ($(size(decoded_bits)) ≂̸ $(size(encoded_bits))"
    last_bit = 0
    transition_bit = (transitions == :high) ? 1 : 0
    for n ∈ eachindex(encoded_bits)
        current_bit = encoded_bits[n]
        if current_bit != last_bit
            decoded_bit = transition_bit
        else
            decoded_bit = 1 - transition_bit
        end
        last_bit = current_bit
        decoded_bits[n] = decoded_bit
    end
    return nothing
end

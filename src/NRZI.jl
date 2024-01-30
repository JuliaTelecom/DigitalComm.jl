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
    encoded_bits = eltype(bits)[]

    last_bit = 0

    transition_bit = (transitions == :high) ? 1 : 0

    for bit in bits
        if bit == transition_bit
            last_bit =  last_bit ⊻ 1
        end
        push!(encoded_bits, last_bit)
    end

    return encoded_bits
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
    decoded_bits = eltype(encoded_bits)[]
    last_bit = 0

    transition_bit = (transitions == :high) ? 1 : 0

    for current_bit in encoded_bits
        if current_bit != last_bit
            decoded_bit = transition_bit
        else
            decoded_bit = 1 - transition_bit
        end
        push!(decoded_bits, decoded_bit)
        last_bit = current_bit
    end

    return decoded_bits
end
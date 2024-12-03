
#= Debugging purpose
println(ZeroPadding(MsgPlusCRC_bits_HammingEncoded, CR_payload, DE_payload))
ZeroPadding( header_bits_HammingEncoded, CR_header, DE_header)
calculate_padding_size_in_bytes(header_bits_HammingEncoded, CR_header, DE_header)
=#

"""
    zero_padding(original_vector::Vector{Int64}, CR::Int8, DE::Bool, SF::Int8) -> Vector{Int64}

Pads a vector of integers with zeros based on the specified coding rate (CR), low data rate optimization (DE), and spreading factor (SF).

# Arguments
- `original_vector`: A vector of integers representing the original data.
- `CR`: Coding rate, an integer (e.g., 1 to 4) defining the redundancy of the Hamming code.
- `DE`: A boolean flag indicating whether Low Data Rate Optimization (LDRO) is enabled.
- `SF`: Spreading Factor, an integer that defines the symbol duration in LoRa modulation.

# Returns
- `Vector{Int64}`: A new vector with zero-padding appended to the original vector.
"""

function zero_padding(original_vector::Vector{Int64}, CR::Int8, DE::Bool, SF::Int8)::Vector{Int64}

    num_of_Bytes = calculate_padding_size_in_bytes(original_vector, CR::Int8, DE::Bool, SF::Int8 )

    # Create a vector of zeros with the specified length
    zeros_vector = zeros(Int64, num_of_Bytes*8)
    
    # Concatenate the original vector with the zeros vector
    padded_vector = vcat(original_vector, zeros_vector)
    return padded_vector
end



#=
vector = MsgPlusCRC_bits_HammingEncoded
original_vector = vector
CR = 3
DE = false
SF
=#

function calculate_padding_size_in_bytes(vector::Vector{Int64}, CR::Int8, DE::Bool, SF::Int8)::Int

    if(DE == true)
        SF_LDRO = SF - 2;
    else 
        SF_LDRO = SF;
    end
   
    # number of bits in the vector
    bit_l = length(vector)
    # number of bytes in the vector
    byte_l = length(vector) / 8
    # the number of bits in a byte minus the bits of the codeword 
    nonImportantZeros = byte_l * (8 - (CR + 4))
    # the codeword bits
    importantBits = bit_l - nonImportantZeros
    # the size of the interleaving block
    Final_blockSize = (CR + 4)  * (SF_LDRO)
    blockSize = Final_blockSize - (importantBits % Final_blockSize)
    
    if (blockSize == Final_blockSize)
        num_of_Bytes = 0
    else
        num_of_Bytes = Int( blockSize  / (CR + 4))
    end
    #blockSize  = Final_blockSize - importantBits
    return num_of_Bytes
end
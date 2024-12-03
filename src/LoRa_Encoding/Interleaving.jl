using Base.Iterators: take, cycle
include("../Utilities.jl")


#=
IntV = fff
CR = 3
DE = false
isheader = false
=#

"""
    interleaving(IntV::Vector{Int64}, CR::Int8, DE::Bool, isheader::Bool, SF::Int8)

This function performs interleaving on an input vector of integers, converting it into a bit matrix and rearranging the bits 
according to specified parameters. It includes support for header interleaving and Low Data Rate Optimization (LDRO).

# Arguments
- `IntV::Vector{Int64}`: The input vector of integers.
- `CR::Int8`: The coding rate parameter.
- `DE::Bool`: Specifies whether LDRO (Low Data Rate Optimization) is enabled.
- `isheader::Bool`: Indicates if the current interleaving process is for a LoRa header.
- `SF::Int8`: Spreading factor used for the interleaving process.

# Returns
- A `BitMatrix` representing the interleaved data.
"""


function interleaving(IntV::Vector{Int64}, CR::Int8, DE::Bool, isheader::Bool, SF::Int8)

    # Convert the input vector of integers to a BitVector
    bV = int64vector_to_bitvector(IntV)

    # Conversion to BitMatrix with 8 rows
    BitM =  bitvector_to_bitmatrix(bV, Int8(8))
    #println(BitM)

    # Adjust spreading factor for Low Data Rate Optimization (LDRO) if enabled
    if(DE == true)
        SF_LDRO = SF - 2;
    else 
        SF_LDRO = SF;
    end

    # Initialize an interleaved matrix of size (CR + 4) x SF_LDRO
    interleaved_matrix = BitMatrix(falses(CR + 4, SF_LDRO)) 

    # Perform interleaving: reverse and circularly shift bits for each row
    for x in 1:(CR + 4)
        println(circshift(reverse(BitM[:, (8 - (CR + 4)) + x]), x-1))
        interleaved_matrix[x, :] = circshift(reverse(BitM[:, (8 - (CR + 4)) + x]), x-1)
    end


     # Extend the matrix for LDRO by appending empty columns if LDRO is enabled
    if(DE == true)
        extended_bit_matrix = BitMatrix(falses(CR + 4, SF))
        # copy matrix into the first SF_LDRO columns
        extended_bit_matrix[:, 1:SF_LDRO] = interleaved_matrix
        interleaved_matrix = extended_bit_matrix
    end

    # If processing a header, calculate and store parity bits for each row
    if(isheader == true)
        for i in 1:size(interleaved_matrix, 1)
            row = interleaved_matrix[i, :]  # Get the i-th row
           # parity_bit = sum(row) % 2  # Calculate even parity (modulo 2)
           # interleaved_matrix[i, end-1] = parity_bit  # Store the parity bit in the 9th column
        end
    end

    return interleaved_matrix
end




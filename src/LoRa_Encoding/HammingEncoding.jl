include("C://Users//User//Desktop//Masterarbeit//Simulation//src//Utilities.jl")



#= Debugging purpose 
bits = BitVector([0,0,1,1,1,1,1,0])
CR_value = Int8(1)

nibbles_hammingEncoded = hamming_encoding(bits,CR_value)

println(nibbles_hammingEncoded)

nibbles_hammingDecoded = hammingDecoding(bv3)
=#

"""
    hamming_encoding(bits::BitVector, CR_value::Int8) -> Vector{Int64}

Encodes a `BitVector` using Hamming encoding based on the specified coding rate (CR). 
Supports CR values of 1 (4/5), 2 (4/6), 3 (4/7), and 4 (4/8).

# Arguments
- `bits`: A `BitVector` representing the data to encode.
- `CR_value`: An integer (1, 2, 3, or 4) indicating the desired coding rate.

# Returns
- `Vector{Int64}`: The Hamming-encoded bit sequence as a vector of integers '0' or '1'
"""


function hamming_encoding(bits::BitVector, CR_value::Int8)::Vector{Int64}
    
        nibbles_hammingEncoded = Int[]

        # cuts the bitvector into pieces of 4 (nibbles)
        nibbles = bitvector_to_pieces(bits, 4)

        # Generator matrix for CR = 2,3,4
        G = BitMatrix([1 1 1 0 1 0 0 0;
                       1 0 1 1 0 1 0 0;
                       0 1 1 1 0 0 1 0;
                       1 1 0 1 0 0 0 1]);

        # Generator matrix for CR = 1
        G2 = BitMatrix([1 1 0 0 0;
                        1 0 1 0 0;
                        1 0 0 1 0 ;
                        1 0 0 0 1]);

        # Iterate over each nibble and encode it based on the CR value
        for i in eachindex(nibbles)
    
            codeword = Int[]  # Initialize an empty container for the codeword

            if (CR_value == 1 ) # CR 4/5
                codeword = mod.(nibbles[i]' * G2, 2); # Multiply nibble by G2 and take modulo 2
                codeword = vec(hcat(0,0,0, codeword)); # Prepend 3 zeros for padding

            elseif (CR_value == 2) # CR 4/6
                codeword = mod.(nibbles[i]' * G[:,3:8], 2) # Multiply nibble by a subset of G
                codeword = vec(hcat(0,0, codeword)); # Prepend 2 zeros for padding

            elseif (CR_value == 3) # CR 4/7
                codeword = mod.(nibbles[i]' * G[:,2:8], 2) # Multiply nibble by subset of G
                codeword = vec(hcat(0, codeword)); # Prepend 1 zero for padding

            elseif (CR_value == 4) # CR 4/8
                codeword = vec(mod.(nibbles[i]' * G[:,1:8], 2)) # Multiply nibble by the full G matrix

            else
                println("Wrong Coding rate @HammingEncoding.jl")
            end

            nibbles_hammingEncoded = [nibbles_hammingEncoded; codeword]

        end 

    return nibbles_hammingEncoded

end
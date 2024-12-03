include("C://Users//User//Desktop//Masterarbeit//Simulation//src//Utilities.jl")

#------- for debugging --------

symbol = [1, 5, 13, 9, 29, 25, 17, 21, 61, 57, 49, 53, 33, 37, 45, 41, 125, 121, 113, 117, 97, 101, 109, 105, 65, 69, 77, 73, 93, 89, 81, 85]

symbol_payload = [196, 104, 499, 498,  393, 191, 163]
symbol_payload_with_error = [196, 105, 499, 498,  393, 193, 163]
symbol_header = [125, 57, 1021, 513, 1, 768, 445, 61]
symbol_header_with_error = [125, 57, 1021, 245, 1, 769, 449, 61]
symbol = symbol_header_with_error
symbol = symbol_payload_with_error

symbol = [105]
DE = false

#xxx = (grayDecoding(symbol, DE, SF))
#=
SF = Int8(10)

xxx_split = bitvector_to_pieces(xxx, 10)
for i in 1:1:length(xxx_split)
    println(xxx_split[i])
    println(bitvector_to_int(xxx_split[i]))
end

println(int_to_bitvector(1021, SF))
println(int_to_bitvector(245, SF))
=#

function grayDecoding(symbol::Vector{Int}, DE::Bool, SF::Int8)::BitVector

    decodedBitVector = BitVector()

    for i in 1:1:length(symbol)

        if (DE == true)
            symbol[i] = map_to_nearest_4k_plus_1(symbol[i])  # Mapping of Symbols for LDRO ON
            symbol[i] += 3        #  add 3 from even elements needed for correct Gray code
            if (mod(symbol[i],8) == 4)      # is symbol not devisable by 8 or not close to it
                symbol[i] -= 3             # substract 3 value
            end
        end

        symbol[i] -= 1

        B_ = BitVector(falses(SF))
        
        B_ = int_to_bitvector(symbol[i], SF) # Chirp Demodulation ich bekomme ein Symbol heraus ein Int Wert
        
        # Gray Coding
        G = BitVector(falses(SF))

        G[1] = B_[1]

        for i in 1:1:length(B_)-1
            G[i+1] = xor(B_[i+1], B_[i])
        end

        decodedBitVector = append!(decodedBitVector, G)
    end

    return decodedBitVector
end

# ------ Mapping of Symbols for LDRO ON ---------

function map_to_nearest_4k_plus_1(value::Int)
    # Berechne das nächstgelegene Symbol der Reihe 4k + 1
    nearest_symbol = 4 * (value ÷ 4) + 1 # ÷ idendical to floor operation   
    return nearest_symbol
end

# Beispielnutzung

array_0_to_127 = collect(0:127)
println(array_0_to_127)
map_to_nearest_4k_plus_1(443)
random_value = 0 # Zufälliger Integer-Wert zwischen 0 und 127
for i in 1:1:128
    println(map_to_nearest_4k_plus_1(array_0_to_127[i]))
end


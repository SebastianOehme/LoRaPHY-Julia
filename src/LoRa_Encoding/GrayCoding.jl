include("C://Users//User//Desktop//Masterarbeit//Simulation//src//Utilities.jl")

# -----------------------------------------------------
# Reverse Gray coding

# the first step is to do a calssical gray decoding to encode the data
# so that one gets the correct integer value for the chirp modulation
# Gray (de)coding

# Anhand dieses Videos erklären https://www.youtube.com/watch?v=MsMkiTcc_w0&t=73s


#println(grayCoding(bitM, true))

bitM = BitMatrix([
    0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 1 0 0;
    0 0 0 0 0 1 0 0 0;
    0 0 0 0 0 1 1 0 0;
    0 0 0 0 1 0 0 0 0;
    0 0 0 0 1 0 1 0 0;
    0 0 0 0 1 1 0 0 0;
    0 0 0 0 1 1 1 0 0;
    0 0 0 1 0 0 0 0 0;
    0 0 0 1 0 0 1 0 0;
    0 0 0 1 0 1 0 0 0;
    0 0 0 1 0 1 1 0 0;
    0 0 0 1 1 0 0 0 0;
    0 0 0 1 1 0 1 0 0;
    0 0 0 1 1 1 0 0 0;
    0 0 0 1 1 1 1 0 0;
    0 0 1 0 0 0 0 0 0;
    0 0 1 0 0 0 1 0 0;
    0 0 1 0 0 1 0 0 0;
    0 0 1 0 0 1 1 0 0;
    0 0 1 0 1 0 0 0 0;
    0 0 1 0 1 0 1 0 0;
    0 0 1 0 1 1 0 0 0;
    0 0 1 0 1 1 1 0 0;
    0 0 1 1 0 0 0 0 0;
    0 0 1 1 0 0 1 0 0;
    0 0 1 1 0 1 0 0 0;
    0 0 1 1 0 1 1 0 0;
    0 0 1 1 1 0 0 0 0;
    0 0 1 1 1 0 1 0 0;
    0 0 1 1 1 1 0 0 0;
    0 0 1 1 1 1 1 0 0;
])

bitM = BitMatrix([
    0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 1 1 0;
    0 0 0 0 0 1 0 1 0;
    0 0 0 0 0 1 1 0 0;
    0 0 0 0 1 0 0 1 0;
    0 0 0 0 1 0 1 0 0;
    0 0 0 0 1 1 0 0 0;
    0 0 0 0 1 1 1 1 0;
    0 0 0 1 0 0 0 1 0;
    0 0 0 1 0 0 1 0 0;
    0 0 0 1 0 1 0 0 0;
    0 0 0 1 0 1 1 1 0;
    0 0 0 1 1 0 0 0 0;
    0 0 0 1 1 0 1 1 0;
    0 0 0 1 1 1 0 1 0;
    0 0 0 1 1 1 1 0 0;
    0 0 1 0 0 0 0 1 0;
    0 0 1 0 0 0 1 0 0;
    0 0 1 0 0 1 0 0 0;
    0 0 1 0 0 1 1 1 0;
    0 0 1 0 1 0 0 0 0;
])


bitM = BitMatrix([
    0 0 0 1 0 0 0 0 0 0;
    0 0 0 0 1 0 0 1 0 0;
    1 0 0 0 0 0 0 0 0 0;
    1 1 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0;
    1 0 1 0 0 0 0 0 0 0;
    0 1 0 1 1 0 0 0 0 0;
    0 0 0 0 1 0 0 0 0 0
])

bitM = BitMatrix([
    0 1 0 1 0 0 0; 
    1 1 1 0 1 0 0;
    1 1 0 1 1 0 0;
    1 1 1 0 1 0 0;
    0 0 0 1 1 0 0
])

bitM = BitMatrix([
    0 0 0 0 1 0 0 1 0 0;
])

# --------- Paper ---------
bitM = BitMatrix([
    0 1 0 1 0 0 0;
    1 1 1 0 1 0 0; 
    1 1 0 1 1 0 0;
    1 1 1 0 1 0 0;
    0 0 0 1 1 0 0;
    0 0 0 1 0 0 0;
    0 0 0 0 0 0 0;
    0 0 1 0 0 0 0
])

bitM = BitMatrix([
    0 1 1 1 0 0;
    0 1 0 1 0 0; 
])

DE = true

#println(grayCoding(bitM, DE))

function grayCoding(bitM::BitMatrix, DE::Bool)::Vector{Int}
      
    row, column = size(bitM)
    #symbol = Vector{Int}(undef, row)
    symbol = Vector{Int}(zeros(row))

    # loop that maps the transmission data to the corresponding integer value for the chirp
    for i in 1:1:length(symbol)
    
        B_ = BitVector(falses(column)) #  stores the gray decoded bitvector 

        # Gray Decoding MSB is the same for the G Gray Code as for the Binary value
        B_[1] = bitM[i,1]
 
        # loop that performs gray decoding
        for k in 1:1:column-1
            B_[k+1] = xor(bitM[i, k+1], B_[k])
        end

        println(B_)
        # Symbol was ich in die Chirp Modulation gebe
        symbol[i] = bitvector_to_int(B_)

       
        # add one to mitigate the bin drift
        symbol[i] += 1

         # LDRO ON
         if (DE == true)
            if (symbol[i] % 2 == 0 && symbol[i] > 2)  # Check if the element is even
                symbol[i] -= 3       # Subtract 3 from even elements
                println("wwwww")
            end 
        end
    end

    # LDRO ON
    #if (DE == true)
    #    for i in 1:length(symbol)
    #        if (symbol[i] % 2 == 0 && symbol[i] > 2)  # Check if the element is even
    #           symbol[i] -= 3       # Subtract 3 from even elements
    #        end
    #    end  
    #end

    return symbol
end



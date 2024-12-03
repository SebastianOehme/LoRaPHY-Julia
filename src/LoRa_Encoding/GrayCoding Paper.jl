

function grayCoding(bitM::BitMatrix, DE::Bool)::Vector{Int}
      
    row, column = size(bitM)
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

        # Symbol was ich in die Chirp Modulation gebe
        symbol[i] = bitvector_to_int(B_)

        # add one to mitigate the bin drift
        symbol[i] += 1
    end

    # LDRO ON
    if (DE == true)
        for i in 1:length(symbol)
            if (symbol[i] % 2 == 0 && symbol[i] > 2)  # Check if the element is even
               symbol[i] -= 3        # Subtract 3 from even elements
            end
        end  
    end

    return symbol
end


# Function to perform Gray reverse coding on a BitMatrix and returns integer symbols
# @param bitM::BitMatrix - Input matrix representing binary-coded data
# @param DE::Bool - Flag indicating if LDRO is enabled
# @return Vector{Int} - Vector of decoded integer symbols after Gray reverse coding and adjustments
function grayReverseCoding(bitM::BitMatrix, DE::Bool)::Vector{Int}
    row, column = size(bitM)               # Get matrix dimensions
    symbol = Vector{Int}(undef, row)       # Initialize vector to store symbols

    # Loop through each row of bitM to calculate each symbol
    for i in 1:row
        # Initialize BitVector to store the Gray decoded bits for each row
        gray_decoded = BitVector(falses(column))

        # Set the Most Significant Bit (MSB) directly (first bit remains the same)
        gray_decoded[1] = bitM[i, 1]

        # Perform Gray decoding for each remaining bit
        for k in 1:column-1
            gray_decoded[k+1] = xor(bitM[i, k+1], gray_decoded[k])
        end

        # Convert the Gray decoded BitVector to an integer and store it in the symbol array
        symbol[i] = bitvector_to_int(gray_decoded) + 1   # Offset to correct binary drift
    end

    # If LDRO is enabled, adjust symbols to handle even values 
    if DE
        for i in 1:length(symbol)
            if symbol[i] % 2 == 0 && symbol[i] > 2  # Check for even symbols greater than 2
               symbol[i] -= 3                       # Adjust by subtracting 3
            end
        end  
    end

    return symbol
end



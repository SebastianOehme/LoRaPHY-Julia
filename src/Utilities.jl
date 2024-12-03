# keine Ahnung ob ich das überhaubt brauche
#using BitArrays 
using LinearAlgebra

function matrix_to_bitvector(bit_matrix)
    bv = vec(bit_matrix')
    return bv
end


function bitvector_to_bitmatrix(bv::BitVector, columns::Int8)
    # Calculate the number of rows
    rows = Int(ceil(length(bv) / columns))
    
    # Pad the BitVector if necessary to ensure the matrix is full
    padded_bv = vcat(bv, falses(rows * columns - length(bv)))
    
    # Reshape the padded BitVector into a BitMatrix with the specified number of columns
    columns = Int64(columns)
    bit_matrix = reshape(padded_bv, columns, rows)
    transposed_matrix = permutedims(bit_matrix)

    return transposed_matrix
end


function bitvector_to_hex(bits::BitVector)::Vector{UInt8}
    n_bytes = div(length(bits), 8)  # Number of bytes to form
    hex_values = Vector{UInt8}(undef, n_bytes)

    for i in 1:n_bytes
        byte = UInt8(0)
        for j in 0:7
            byte |= (bits[8 * (i - 1) + (8 - j)] ? 1 : 0) << j
        end
        hex_values[i] = byte
    end
    
    return hex_values
end





"""
    hex_to_bitvector(hex_values::Vector{UInt8}) -> BitVector

Converts a vector of hexadecimal (UInt8) values into a `BitVector` representing the corresponding binary bits.

# Arguments
- `hex_values`: A vector of `UInt8` values to convert to a `BitVector`.

# Returns
- A `BitVector` where each byte in `hex_values` is expanded into its 8-bit binary representation.

# Example
```julia
hex_values = [0xA3, 0x1F]
bits = hex_to_bitvector(hex_values)
println(bits)  # Output: BitVector[1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1]
"""


function hex_to_bitvector(hex_values::Vector{UInt8})::BitVector 
    bits = BitVector(undef, 8 * length(hex_values))
    for (i, byte) in enumerate(hex_values)
        for j in 0:7
            bits[8 * (i - 1) + (8 - j)] = (byte >> j) & 1 == 1
        end
    end
    return bits
end




# möglicherweise sollte der Int64 vektor nur aus 0 und 1 bestehen damit es funktioniert
function int64vector_to_bitvector(vec::Vector{Int64})::BitVector
    bitvec = BitVector()
    for num in vec
        push!(bitvec, num == 1)
    end
    return bitvec
end


function bitvector_to_pieces(bv::BitVector, n::Int)
    pieces = Vector{BitVector}(undef, 0)  # Initialize an empty array to store BitVectors
    for i in 1:n:length(bv)
        piece = BitVector(undef, n)  # Create a BitVector of length 4
        for j in 0:(n-1)
            if i + j <= length(bv)
                piece[j + 1] = bv[i + j]
            else
                piece[j + 1] = false  # Fill remaining bits with false (0) if not enough bits
            end
        end
        pieces = push!(pieces, piece)
    end
    return pieces
end

function bitvector_to_int(bits::BitVector)::Int
    # Calculate the integer value directly from the BitVector
    return sum(b << (length(bits) - i - 1) for (i, b) in enumerate(bits))
end

function bitvector_to_int(bits::BitVector)::Int
    # Calculate the integer value directly from the BitVector
    return sum(bits[i] * 2^(length(bits) - i) for i in 1:length(bits))
end

function int_to_bitvector(x::Int, length::Int8)::BitVector
    bits = falses(length)  # Create a BitVector of the specified length initialized to false (0)
    for i in 1:length
        bits[length - i + 1] = (x % 2 == 1)  # Check the least significant bit
        x = x >> 1  # Right shift the integer to check the next bit
    end
    return bits
end

# ------------- Not used in main programm  -------------------


function find_rise_fall_indices(arr)
    n = length(arr)
    start_idx = 0
    end_idx = 0
    
    # Find where the sequence starts (rising after a non-positive value)
    for i in 2:n
        if arr[i-1] <= 0 && arr[i] > 0
            start_idx = i
            break
        end
    end
    
    # Check if a valid start index was found
    if start_idx == 0
        return nothing  # No valid sequence found
    end
    
    # Find where the sequence ends (falling but last positive element included)
    for i in start_idx+1:n
        if arr[i] <= 0 || (arr[i] < arr[i-1] && arr[i] <= 0)
            end_idx = i - 1
            break
        end
        # Include the last positive value even if it's falling
        if i == n
            end_idx = i
        end
    end
    
    # Handle case where the array ends with a positive value
    if end_idx == 0
        end_idx = n
    end
    
    return start_idx,end_idx
end

# Example array
#arr = [100, 44, 32, -9, -12, -5, 1, 6, 89, 55, 2, -6, -8, 2]

# Find the indices of the rising and falling sequence
#i,j = find_rise_fall_indices(arr)


# Original array
#array = magnitude_of_symbols
#window_size 
#index = S_demodulated

# Function to calculate the padded window
function padded_window(array, window_size, index)
    left_pad = 0
    right_pad = 0
    half_window = (window_size - 1) ÷ 2

    # Start and end indices of the window
    start_index = index - half_window
    end_index = index + half_window

    # If start_index is less than 1, we need to pad on the left
    if start_index < 1
        left_pad = abs(start_index)+1
        # Get the left padding by cycling the elements from the end
        padding = array[end-left_pad+1:end]
        padded_array = vcat(padding, array)
      
    elseif end_index > length(array)
        # If end_index exceeds the length, pad on the right side
        right_pad = end_index - length(array)
        # Get the right padding by cycling the elements from the beginning
        padding = array[1:right_pad]
        padded_array = vcat(array, padding)
    else
        padded_array = array
    end

    # Extract the window from the padded array
    window = padded_array;
    
    return window, left_pad, right_pad
end


# Get the padded window
#window,l,r = padded_window(array, window_size, index)





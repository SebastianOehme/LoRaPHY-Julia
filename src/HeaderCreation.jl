"""
    header_creation(msg_bytes::Vector{UInt8}, CRC_enable::Bool, CR_payload::Int8) -> BitVector

Creates a LoRa header in the form of a `BitVector`, encoding the payload length, coding rate (CR), CRC flag, 
and the calculated CRC for the header itself.

# Arguments
- `msg_bytes`: A vector of `UInt8` values representing the message payload.
- `CRC_enable`: A boolean indicating whether the CRC is enabled for the payload.
- `CR_payload`: The coding rate of the payload, specified as an integer (1,2,3 or 4).

# Returns
- `BitVector`: The constructed header as a `BitVector` with a total length of 20 bits.
"""

function header_creation(msg_bytes::Vector{UInt8}, CRC_enable::Bool, CR_payload::Int8)::BitVector

    # Initialize the header
    header = BitVector(undef,20)
    payload_length = UInt8(length(msg_bytes))

    # Store payload_length into positions 1 to 8 of the BitVector
    set_payload_bits!(header, payload_length, 1)

    # Store the CR bits into positions 9 to 11 of the BitVector
    set_cr_bits!(header, CR_payload , 9)

    # Set CRC_enable
    header[12] = CRC_enable;

    # Calculate Header CRC bits
    crc_header_bits = calc_crc_header(header[1:12]);

    # Set Header CRC bits
    header[16:20] = crc_header_bits;

    return header;
end




# calculates the header CRC
function calc_crc_header(header::BitVector)::BitVector

    h = header'
    M = BitMatrix([ 1 1 0 0 0;
                    1 0 1 0 0;
                    1 0 0 1 0;
                    1 0 0 0 1;
                    0 1 1 0 0;
                    0 1 0 1 0;
                    0 1 0 0 1;
                    0 0 1 1 0;
                    0 0 1 0 1;
                    0 0 0 1 1;
                    0 0 1 1 1;
                    0 1 0 1 1])

    crc_vector = mod.(h*M,2)

    #converts the crc_vector from a Int64 Matrix to a BitVector
    return BitVector(vec(crc_vector))

end


# Function to set bits in BitVector from a UInt8 value
function set_payload_bits!(header::BitVector, value::UInt8, start_pos::Int)
    for i in 0:7
        header[start_pos + i] = (value >> (7 - i)) & 1 == 1
    end
end

# Function to set bits in BitVector from a 3-bit value
function set_cr_bits!(header::BitVector, value::Int8, start_pos::Int)
    for i in 0:2
        if start_pos + i <= length(header)  # Prevent out-of-bounds errors
            header[start_pos + i] = (value >> (2 - i)) & 1 == 1
        end
    end
end









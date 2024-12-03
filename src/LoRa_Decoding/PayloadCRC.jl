include("C://Users//User//Desktop//Masterarbeit//Simulation//src//Utilities.jl")

# Define CRC-16 Polynomial
const CRC16POLY = 0x1021


function checkPayloadwithCRC(msg_bits::BitVector, crc_bits::BitVector)

    # Define the bitstream
    msg_and_crc_bits = vcat(msg_bits, crc_bits);
 
    bitcount = length(msg_and_crc_bits)

    # Initialize the result
    result = UInt16(0)

    # Calculate the CRC-16
    for i in 1:bitcount
        if ((result >> 15) & 1) != msg_and_crc_bits[i]
            result = (result << 1) ⊻ CRC16POLY
        else
            result <<= 1
        end
    end

    if (result == 0)
        print("No Error Detectet from Payload CRC")
    else
        print("Error Detectet from Payload CRC!!!")
    end

end
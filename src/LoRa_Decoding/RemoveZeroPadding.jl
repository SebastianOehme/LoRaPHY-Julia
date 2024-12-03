function removeZeroPadding(bit_M::BitMatrix, isheader::Bool, numberOfNibblesPayload::Int64, CRC_enabled_decoded::Bool)

    if (isheader == true)
        bit_M = bit_M[1:5,:];
    else
        if (CRC_enabled_decoded == true)
            CRC_nibbles = 4
        else
            CRC_nibbles = 0 
        end
        bit_M = bit_M[1:(numberOfNibblesPayload + CRC_nibbles),:];
    end

    return bit_M
end

TOA(51)

#SF = 12
#Bw = 31250
#IH = 0
#CRC_enable = 1
#CR_payload = 4
#DE_payload = 1
#PL = 200

#preamble_length = 8


function TOA(PL::Int64)
    # Symbol Time
    T_s =  (2^SF)/Bw

    # Time for preamble
    T_preamble = (preamble_length + 4.25) * T_s

    # number of payload symbols + header etc. 
    n_payload = 8 + max(ceil( (8* PL - 4*SF + 28 + 16* CRC_enable - 20 * IH) / (4 * (SF - 2 * DE_payload)) ) * (CR_payload + 4), 0)

    # Time for payload symbols + header etc. 
    T_payload = n_payload * T_s

    # total time on air of the packet in [s]
    T_packet = T_preamble + T_payload

    return T_packet
end
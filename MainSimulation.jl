# List of imported Packets
using Plots


include("src/Utilities.jl")

include("src/HeaderCreation.jl")
include("src/LoRa_Encoding/PayloadCRC.jl")

# Encoding
include("src/LoRa_Encoding/Whitening.jl")
include("src/LoRa_Encoding/HammingEncoding.jl")
include("src/LoRa_Encoding/ZeroPadding.jl")
include("src/LoRa_Encoding/Interleaving.jl")
include("src/LoRa_Encoding/GrayCoding.jl")

# Decoding
include("src/LoRa_Decoding/GrayDecoding.jl")
include("src/LoRa_Decoding/Deinterleaving.jl")
include("src/LoRa_Decoding/RemoveZeroPadding.jl")
include("src/LoRa_Decoding/HammingDecoding.jl")
include("src/LoRa_Decoding/Dewhitening.jl")
include("src/LoRa_Decoding/HeaderCRC.jl")

#include("src/LoRa_Decoding/PayloadCRC.jl")


# Modulation/Demodulation
#include("src/LoRa_Modulation/LoRaModulation.jl")
#include("src/LoRa_Modulation/ChirpSignal.jl")
#include("src/LoRa_Demodulation/Dechirping.jl")


# TODO 
# 1. Implicit header einbauen
# 2. CRC_enabeld = false das dann auch CRC nicht eingebaut wird!
# 3. Wenn man die Function LoRa() verwendet sollten zum Beipiel bei createHeader die Werte Übergeben werden
# 4. Testen
# 5. LoRa packete enabelen
# 6. Alle functionen mitRückgabewert()::Ausstatten
# 7. Write comments that explain "why," not "what," unless the code is complex.
# 8. Hamming encoding returns Integer 0 and 1's could be changed to bits if possible

   # main(msg)


#function main(msg::String)


SF = Int8(10);              # {7,8,9,10,11,12}
CR_payload = Int8(3);       # CR = 1 (4/5) ; CR = 2 (4/6) ; CR = 3 (4/7) ; CR = 4 (4/8) 
CRC_enable = true;    # {true, false}
DE_payload = false;   # {true, false} is true when T_s is > 16 ms
IH = 0;               # {0,1} when header is enabled 0 when not 1
# -------------------------------------------
preamble_length = 8;        # {6-65535}
payload_length_fixed = 50;  # {0-255}

# fixed LoRa parameters
CR_header = Int8(4);
DE_header = true; 

# Define the message
msg = "Hi"


# Test setup
DE_payload = false;
CR_payload = Int8(1);       # CR = 1 (4/5) ; CR = 2 (4/6) ; CR = 3 (4/7) ; CR = 4 (4/8) 

SF = Int8(7)
LoRa(msg, SF, CR_payload, DE_payload, CRC_enable, CR_header, DE_header)
SF = Int8(8)
LoRa(msg, SF, CR_payload, DE_payload, CRC_enable, CR_header, DE_header)
SF = Int8(9)
LoRa(msg, SF, CR_payload, DE_payload, CRC_enable, CR_header, DE_header)
SF = Int8(10)
LoRa(msg, SF, CR_payload, DE_payload, CRC_enable, CR_header, DE_header)
SF = Int8(11)
LoRa(msg, SF, CR_payload, DE_payload, CRC_enable, CR_header, DE_header)
SF = Int8(12)
LoRa(msg, SF, CR_payload, DE_payload, CRC_enable, CR_header, DE_header)



#function LoRa(msg::String, SF::Int8, CR_payload::Int8, DE_payload::Bool, CRC_enable::Bool, CR_header::Int8,  DE_header::Bool)::String

   
    # --------------------- Message Creation -------------------

  

    # Convert the string to bytes
    msg_bytes = [UInt8(c) for c in msg]

    # bytes to bits
    msg_bits = hex_to_bitvector(msg_bytes)
    println(msg_bits) 
    SF
    # --------------------------- HEADER ENCODING ----------------------------------    

    # --------------------------- Header Creation ----------------------------------
    header_bits = header_creation(msg_bytes, CRC_enable, CR_payload) # BitVector
    println(header_bits)
SF

    # ----------------------- Hamming Encoding - Header -----------------------------    
    header_bits_HammingEncoded = hamming_encoding(header_bits, CR_header)
    bitvector_to_bitmatrix(int64vector_to_bitvector(header_bits_HammingEncoded),8)

SF
    
    # -------------------------- Zero Padding - Header ------------------------------
    header_bits_zeroPadded = zero_padding(header_bits_HammingEncoded, CR_header, DE_header, SF)
    for i in 1:8:length(header_bits_zeroPadded)-7
        println(header_bits_zeroPadded[i:i+7])
    end
SF
    
    # --------------------------- Interleaving - Header -----------------------------
    header_interleaved_Matrix = interleaving(header_bits_zeroPadded, CR_header, DE_header, true, SF)
    println(header_interleaved_Matrix)

    # ----------------------- Gray "reverse" coding - Header  -----------------------
    symbols_header = grayCoding(header_interleaved_Matrix, DE_header)
    println(symbols_header)
        


    # ------------------------------- PAYLOAD ENCODING -------------------------------    

    # ---------------------------- CRC Calculation - Payload -------------------------
    crc_bytes = PacketComputeCrc(buffer, CRC_TYPE_CCITT)
    println(hex_to_bitvector(crc_bytes))
SF

    # ----------------------------- Whitening - Payload ------------------------------  
    msg_bytes_whithened = whitening(msg_bytes)
    print(hex_to_bitvector(msg_bytes_whithened))
SF

    # ----------------------------- Hamming - Payload + CRC --------------------------
    MsgPlusCRC_bytes = vcat(msg_bytes_whithened, crc_bytes)
    MsgPlusCRC_bits = hex_to_bitvector(MsgPlusCRC_bytes)
    println(MsgPlusCRC_bits)

    MsgPlusCRC_bits_HammingEncoded = hammingEncoding(MsgPlusCRC_bits, CR_payload)
    for i in 1:8:length(MsgPlusCRC_bits_HammingEncoded)-7
        println(MsgPlusCRC_bits_HammingEncoded[i:i+7])
    end

    
    SF
    # ------------------------------ Zero Padding - Payload -------------------------- 
    MsgPlusCRC_bits_HammingEncoded_bV = int64vector_to_bitvector(MsgPlusCRC_bits_HammingEncoded)

    MsgPlusCRC_bits_zeroPadded = ZeroPadding(MsgPlusCRC_bits_HammingEncoded, CR_payload, DE_payload, SF)

    for i in 1:8:length(MsgPlusCRC_bits_zeroPadded)-7
        println(MsgPlusCRC_bits_zeroPadded[i:i+7])
    end


    # ------------------ For Loop necessary for longer Messages ----------------------- 
    length(MsgPlusCRC_bits_zeroPadded)
    # if message is longer than (CR + 4)  * (SF_LDRO) for loop

    if(DE_payload == true)
        SF_LDRO = SF - 2;
    else 
        SF_LDRO = SF;
    end

    blockSize = 8  * (SF_LDRO)
    #blockSize = (CR_payload + 4) * SF_LDRO
    symbols_payload = Int64[]
    length(MsgPlusCRC_bits_zeroPadded)
    SF

    for i in 0:blockSize:length(MsgPlusCRC_bits_zeroPadded)-1
        #i = 0
        #fff = MsgPlusCRC_bits_zeroPadded[i+1:i+blockSize]
        # ---------------------------- Interleaving - Payload + CRC ----------------------
        payload_interleaved_Matrix = interleaving(MsgPlusCRC_bits_zeroPadded[i+1:i+blockSize], CR_payload, DE_payload, false, SF);
        #println(payload_interleaved_Matrix)
        
        # ------------------------- Gray "reverse" coding - Payload ----------------------
        symbols_payload_block = grayCoding(payload_interleaved_Matrix, DE_payload);
        symbols_payload = vcat(symbols_payload, symbols_payload_block)
    end

    #payload_interleaved_Matrix
    #symbols_payload

    # --------------- concatinating Header symbols with payload symbols --------------- 
    symbols = vcat(symbols_header, symbols_payload)
    println(symbols)
    symbols .+= 1
    println(symbols)



    # TODO Je nach Payload length muss man hier cutten  -> will ich nicht machen, packet wird dann halt lang
    preamble = zeros(Int, preamble_length)
 
    LoRaPacket = vcat(preamble, symbols)

    # Der orbit wurde vom user oben definiert. 
    # 





    # --------------------------------------------------------------------
    # Vollständiges Signal mit Orbitsimulation paaren
    # die Orbitsimulation gibt den dopplershift zu ganz genauen Zeitpunkten zurück 
    # somit kann man sich dann den Einfluss auf ein LoRa packet anschauen,
    # im Verlauf eines Orbits 
    # --------------------------------------------------------------------

    
    #S = dechirpSignal(LoRaPacket)
    #print(S)


    # --------------------------------------------------------------------
    # ------------------------- DECODING ---------------------------------

    # TODO Decode packet per Packet


    # --------------------------- HEADER DECODING ----------------------------------    
    
    symbols[1:8]
    # ------------------------ Gray Decoding - Header -------------------------------   
    # Number of Header Symbols is always 8
    header_bV = grayDecoding(symbols[1:8], DE_header, SF)
    header_M = bitvector_to_bitmatrix(header_bV, SF)

    bitvector_to_int(header_M[1,:])
    bitvector_to_int(header_M[2,:])
    bitvector_to_int(header_M[3,:])
    bitvector_to_int(header_M[4,:])
    bitvector_to_int(header_M[5,:])
    bitvector_to_int(header_M[6,:])
    bitvector_to_int(header_M[7,:])
    bitvector_to_int(header_M[8,:])

    SF
    # ----------------------- Deinterleaving - Header -------------------------------   
    deinterleaved_matrix_header = diag_deinterleaving(header_M, CR_header, DE_header)
    deinterleaved_matrix_header
    SF
    # ----------------------- Remove Padding - Header -------------------------------   
    deinterleaved_matrix_header = removeZeroPadding(deinterleaved_matrix_header, true, 0, false)
    SF
    # ----------------------- Hamming Decoding - Header -----------------------------    
    hammingdecoded_header = hammingDecoding(deinterleaved_matrix_header, CR_header)
    println(hammingdecoded_header)
    SF
    # -------------------------- CRC Check header ------------------------------
    check_header_with_crc(hammingdecoded_header)

    SF
    # --------------------- Decode Header Information --------------------------
    payload_length_bytes_decoded = bitvector_to_int(hammingdecoded_header[1:8])
    CR_Payload_decoded = bitvector_to_int(hammingdecoded_header[9:11])
    CRC_enabled_decoded = hammingdecoded_header[12]

    # Print out the Information
    print("Payload length in bytes: " , bitvector_to_int(hammingdecoded_header[1:8]))
    print("Coding Rate: 4/" , bitvector_to_int(hammingdecoded_header[9:11]) + 4)
    print("CRC Enabled: " , hammingdecoded_header[12])

SF

    # ---------------------------- PAYLOAD DECODING ---------------------------------    
    
    symbols[num_header_symbols+1:end]

    # ------------------------ Gray Decoding - Payload -------------------------------   
    payload_bV = grayDecoding(symbols[num_header_symbols+1:end], DE_payload, SF)
    payload_M = bitvector_to_bitmatrix(payload_bV, SF)
    println(payload_M)
    bitvector_to_int(payload_M[1,:])
    bitvector_to_int(payload_M[2,:])
    bitvector_to_int(payload_M[3,:])
    bitvector_to_int(payload_M[4,:])
    bitvector_to_int(payload_M[5,:])
    bitvector_to_int(payload_M[6,:])
    bitvector_to_int(payload_M[7,:])





    
    SF
    # workaround the gray decode problem -3 + 1
    if (DE_payload == true)
        payload_M[:, SF - 1] .= 0
    end

    # ------------------ For Loop necessary for longer Messages ----------------------- 
   
    
    
    # ----------------------- Deinterleaving - Payload -------------------------------   
    row_n,column_n = size(payload_M)
    number_of_loop_passages = Int(row_n / (CR_payload+4))
    deinterleaved_matrix_payload_full =  BitMatrix(undef, 0, 8)
    SF
    for i in 0:1:number_of_loop_passages-1
        payload_M_block = payload_M[(i*(CR_payload+4)+1):(CR_payload+4)*(i+1),:]

        deinterleaved_matrix_payload = diag_deinterleaving(payload_M_block, CR_payload, DE_payload)

        deinterleaved_matrix_payload_full =  vcat(deinterleaved_matrix_payload_full, deinterleaved_matrix_payload)
    end
    SF
    deinterleaved_matrix_payload_full
    SF
    # ----------------------- Remove Padding - Payload ------------------------------- 
    NumberOfPayloadNibbles = payload_length_bytes_decoded*2
    deinterleaved_matrix_payload_full = removeZeroPadding(deinterleaved_matrix_payload_full, false, NumberOfPayloadNibbles , CRC_enabled_decoded)

    CR_payload

    # ----------------------- Hamming Decoding - Payload -----------------------------    
    hammingdecoded_payload = hammingDecoding(deinterleaved_matrix_payload_full, CR_payload)
    println(hammingdecoded_payload)

    # -------------------------- Dewhitening - Payload -------------------------------
    msg_hex = dewhitening(hammingdecoded_payload[1:payload_length_bytes_decoded*8], false)
    println(msg_hex)
    
    # ---------------------------- CRC Check - Payload -------------------------------
    # extract the CRC values
    crc_bits = hammingdecoded_payload[(payload_length_bytes_decoded*8 + 1):end]
    crc_hex = bitvector_to_hex(crc_bits)
    msg_plus_crc_hex = vcat(msg_hex, crc_hex)
    PacketComputeCrc(msg_plus_crc_hex, true)


    # --------------------------- Payload bits to text -------------------------------
    msgDecoded = String(msg_hex)


    # ----------------------------- END -----------------------------------

    return msgDecoded

#end




# --------------- SHORT VERSION OF FUNCTION LORA ------------------------

function LoRa(msg::String, SF::Int8, CR_payload::Int8, DE_payload::Bool, CRC_enable::Bool, CR_header::Int8,  DE_header::Bool)::String

   
    # --------------------- Message Creation -------------------
    # Convert the string to bytes
    msg_bytes = [UInt8(c) for c in msg];
    # bytes to bits
    msg_bits = hex_to_bitvector(msg_bytes);
 
    # --------------------------- HEADER ENCODING ----------------------------------    

    # --------------------------- Header Creation ----------------------------------
    header_bits = headerCreation(msg_bytes, CRC_enable, CR_payload); # BitVector

    # ----------------------- Hamming Encoding - Header -----------------------------    
    header_bits_HammingEncoded = hammingEncoding(header_bits, CR_header);
    
    # -------------------------- Zero Padding - Header ------------------------------
    header_bits_zeroPadded = ZeroPadding(header_bits_HammingEncoded, CR_header, DE_header, SF);
    
    # --------------------------- Interleaving - Header -----------------------------
    header_interleaved_Matrix = interleaving(header_bits_zeroPadded, CR_header, DE_header, true, SF);

    # ----------------------- Gray "reverse" coding - Header  -----------------------
    symbols_header = grayCoding(header_interleaved_Matrix, DE_header);


    # ------------------------------- PAYLOAD ENCODING -------------------------------    

    # ---------------------------- CRC Calculation - Payload -------------------------
    crc_bytes = PacketComputeCrc(msg_bytes, false);

    # ----------------------------- Whitening - Payload ------------------------------  
    msg_bytes_whithened = LoRa_encode_white(msg_bytes, false);

    # ----------------------------- Hamming - Payload + CRC --------------------------
    MsgPlusCRC_bytes = vcat(msg_bytes_whithened, crc_bytes);
    MsgPlusCRC_bits = hex_to_bitvector(MsgPlusCRC_bytes);

    MsgPlusCRC_bits_HammingEncoded = hammingEncoding(MsgPlusCRC_bits, CR_payload);
 
    # ------------------------------ Zero Padding - Payload -------------------------- 
    MsgPlusCRC_bits_zeroPadded = ZeroPadding(MsgPlusCRC_bits_HammingEncoded, CR_payload, DE_payload, SF);

    # ------------------ For Loop necessary for longer Messages ----------------------- 
    if(DE_payload == true)
        SF_LDRO = SF - 2;
    else 
        SF_LDRO = SF;
    end

    blockSize = 8  * (SF_LDRO);
    symbols_payload = Int64[]
   
    for i in 0:blockSize:length(MsgPlusCRC_bits_zeroPadded)-1
       
        # ---------------------------- Interleaving - Payload + CRC ----------------------
        payload_interleaved_Matrix = interleaving(MsgPlusCRC_bits_zeroPadded[i+1:i+blockSize], CR_payload, DE_payload, false, SF);

        # ------------------------- Gray "reverse" coding - Payload ----------------------
        symbols_payload_block = grayCoding(payload_interleaved_Matrix, DE_payload);
        symbols_payload = vcat(symbols_payload, symbols_payload_block);
    end

    # --------------- concatinating Header symbols with payload symbols --------------- 
    symbols = vcat(symbols_header, symbols_payload);
       


    # --------------------------------------------------------------------
    # ------------------------- DECODING ---------------------------------


    # --------------------------- HEADER DECODING ----------------------------------    
    
    # ------------------------ Gray Decoding - Header -------------------------------   
    # Number of Header Symbols is always 8
    header_bV = grayDecoding(symbols[1:8], DE_header, SF);
    header_M = bitvector_to_bitmatrix(header_bV, SF);

    # ----------------------- Deinterleaving - Header -------------------------------   
    deinterleaved_matrix_header = diag_deinterleaving(header_M, CR_header, DE_header);

    # ----------------------- Remove Padding - Header -------------------------------   
    deinterleaved_matrix_header = removeZeroPadding(deinterleaved_matrix_header, true, 0, false);

    # ----------------------- Hamming Decoding - Header -----------------------------    
    hammingdecoded_header = hammingDecoding(deinterleaved_matrix_header, CR_header);
 
    # -------------------------- CRC Check header ------------------------------
    check_header_with_crc(hammingdecoded_header);

    # --------------------- Decode Header Information --------------------------
    payload_length_bytes_decoded = bitvector_to_int(hammingdecoded_header[1:8]);
    CR_Payload_decoded = bitvector_to_int(hammingdecoded_header[9:11]);
    CRC_enabled_decoded = hammingdecoded_header[12];

    # Print out the Information
    println("Payload length in bytes: " , bitvector_to_int(hammingdecoded_header[1:8]))
    println("Coding Rate: 4/" , bitvector_to_int(hammingdecoded_header[9:11]) + 4)
    println("CRC Enabled: " , hammingdecoded_header[12])

    # ---------------------------- PAYLOAD DECODING ---------------------------------    
    
    # ------------------------ Gray Decoding - Payload -------------------------------   
    payload_bV = grayDecoding(symbols[num_header_symbols+1:end], DE_payload, SF);
    payload_M = bitvector_to_bitmatrix(payload_bV, SF);
    
    # workaround the gray decode problem -3 + 1
    if (DE_payload == true)
        payload_M[:, SF - 1] .= 0;
    end

    # ------------------ For Loop necessary for longer Messages ----------------------- 
   
    
    # ----------------------- Deinterleaving - Payload -------------------------------   
    row_n,column_n = size(payload_M);
    number_of_loop_passages = Int(row_n / (CR_payload+4));
    deinterleaved_matrix_payload_full =  BitMatrix(undef, 0, 8);
  
    for i in 0:1:number_of_loop_passages-1
        payload_M_block = payload_M[(i*(CR_payload+4)+1):(CR_payload+4)*(i+1),:];

        deinterleaved_matrix_payload = diag_deinterleaving(payload_M_block, CR_payload, DE_payload);

        deinterleaved_matrix_payload_full =  vcat(deinterleaved_matrix_payload_full, deinterleaved_matrix_payload);
    end
  
    # ----------------------- Remove Padding - Payload ------------------------------- 
    NumberOfPayloadNibbles = payload_length_bytes_decoded*2;
    deinterleaved_matrix_payload_full = removeZeroPadding(deinterleaved_matrix_payload_full, false, NumberOfPayloadNibbles , CRC_enabled_decoded);


    # ----------------------- Hamming Decoding - Payload -----------------------------    
    hammingdecoded_payload = hammingDecoding(deinterleaved_matrix_payload_full, CR_payload);

    # -------------------------- Dewhitening - Payload -------------------------------
    msg_hex = dewhitening(hammingdecoded_payload[1:payload_length_bytes_decoded*8], false);
    
    # ---------------------------- CRC Check - Payload -------------------------------
    # extract the received CRC value
    crc_bits_received = hammingdecoded_payload[(payload_length_bytes_decoded*8 + 1):end];
    crc_hex = bitvector_to_hex(crc_bits_received);

    # calculate the CRC value with the received msg
    crc_bytes_calc_with_received = PacketComputeCrc(msg_hex, true);
    if (crc_hex[1] == crc_bytes_calc_with_received[1] && crc_hex[2] == crc_bytes_calc_with_received[2])
        println("Payload error free")
    else
        println("Error in payload")
    end




    # --------------------------- Payload bits to text -------------------------------
    msgDecoded = String(msg_hex);

    # ----------------------------- END -----------------------------------

    return msgDecoded

end




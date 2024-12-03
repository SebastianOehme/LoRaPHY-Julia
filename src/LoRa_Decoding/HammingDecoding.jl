include("C://Users//User//Desktop//Masterarbeit//Simulation//src//Utilities.jl")

#= for Debugging
bits = BitVector([1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1,
0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1]);

hammingDecoding(bits, 1)


CR_value = 4
bit_M = deinterleaved_matrix_payload_full

row = bit_M[3,:]
println(row)


=#

function hammingDecoding(bit_M::BitMatrix, CR_value::Int8)

    nibbles_hammingDecoded = Vector{BitVector}(undef, 0)
    num_errors = true;

    for row in eachrow(bit_M)

        if (CR_value == 4)
            # 2:8 must match th e 1
            p4 = sum(row[2:8]) % 2  # Calculate even parity (modulo 2)
            if (p4 == row[1])
                # Zero or even Number of errors
                # Continue with CR_Value = 3
                # If this results in an error 2 errors detected!!! Can not be corrected must be send again
                num_errors = false;
                CR_value = 3;
            else 
                # drop p4 and continue with CR_value = 3"
                println("Odd number of error most probably 1 have been detected")
                num_errors = true;
                CR_value = 3;
            end
        end

        if (CR_value == 3)

            # reorder wie aufgeschrieben
            codeword_reordered = BitVector(undef, 7)
            codeword_reordered[4] = row[2]
            codeword_reordered[2] = row[3]
            codeword_reordered[1] = row[4]
            codeword_reordered[6] = row[5]
            codeword_reordered[3] = row[6]
            codeword_reordered[7] = row[7]
            codeword_reordered[5] = row[8]
            

            H = BitMatrix([0 0 1;
                           0 1 0;
                           0 1 1;
                           1 0 0;
                           1 0 1;
                           1 1 0;
                           1 1 1]);

            # error position ausrechnen
            syndrom = BitVector(vec(mod.(codeword_reordered' * H,2)))

            # which is the exact position (gives the position of the error in codeword_reordered)
            position = sum(syndrom[i] * 2^(length(syndrom) - i) for i in 1:length(syndrom))

            # if zero or even number of error where detected with the parity bit p4 which is a checksum 
            # but an error was detected with the Hamming (7,4) than one can conclude that two errors must 
            # have occuerd, in order to prevent that a correct bit is flipped a if else statement is used
            if (num_errors == false && position != 0)
                    println("Two Errors in one nibble has been detected!!!")
                    # TODO here a flag that the packet shoud be send again could be helpful???
            else
                # flip the incorrect bit @ index position 
                if position != 0
                    codeword_reordered[position] = !codeword_reordered[position]
                    println("1 error has been corrected")
                end
            end

            # get rid of the parity bits an reorder LSB right
            nibble = BitVector(undef, 4)
            nibble[1] = codeword_reordered[6]
            nibble[2] = codeword_reordered[3]
            nibble[3] = codeword_reordered[7]
            nibble[4] = codeword_reordered[5]
                    
            nibbles_hammingDecoded = push!(nibbles_hammingDecoded, nibble)
        end

            
        if (CR_value == 2)

            H = BitMatrix([1 0;
                           0 1;
                           1 0;
                           1 1;
                           1 1;
                           0 1]);
    
            syndrom  = mod.(row[3:8]' * H, 2)
    
            if (syndrom[1] != 0 || syndrom[2] != 0)
                println("One Bit error in one Nibble")
            end
    
            nibble = BitVector(undef, 4)
            nibble = row[5:8]
            nibbles_hammingDecoded = push!(nibbles_hammingDecoded, nibble) 
    
        end 

        if (CR_value == 1)

            p5 = sum(row[5:8]) % 2  # Calculate even parity (modulo 2)

            if (p5 != row[4])
                println("One Bit error in on Nibble")
            end

            nibble = BitVector(undef, 4)
            nibble = row[5:8]
            nibbles_hammingDecoded = push!(nibbles_hammingDecoded, nibble)

        end
    end

    return reduce(vcat, nibbles_hammingDecoded)
    
end
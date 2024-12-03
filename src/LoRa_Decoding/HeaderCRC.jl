

function check_header_with_crc(Header_BitVector)

    h = Header_BitVector;

    M =[1 1 0 0 0;
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
        0 1 0 1 1];

    crc_vector = mod.(h[1:12]'*M,2)


    #check if all crc_vector is similar to the transmitted crc values
    bv = crc_vector[1:5] .!= h[16:20];

    # Check if all elements are zeros (false)
    all_zeros = all(x -> x == false, bv)

    if all_zeros == true
        println("no error in header")
    else
        println("ERROR in header!")
    end
end


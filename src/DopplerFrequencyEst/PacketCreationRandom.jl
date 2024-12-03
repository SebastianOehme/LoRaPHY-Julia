function fill_random_packet(packet_size::Int, preamble::Int, SF::Int)

    # Initialize an array of zeros of the given size
    packet = zeros(Int, packet_size)

    # Fill the remaining part of the array with random integers from 0 to 2^SF - 1
    packet[preamble+1:end] .= rand(0:2^SF-1, packet_size - preamble)

    return packet
end



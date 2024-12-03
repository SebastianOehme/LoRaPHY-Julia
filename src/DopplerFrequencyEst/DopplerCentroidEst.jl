fd = 7000 / c * fc


window_size = 40
padding_factor = 4

magnitude_of_symbols = performFFTPadded(dechirpWindowing(chirp_signals[10]), 4*4096)
plot(magnitude_of_symbols)

 #dopplerEstimationCentroidCorrectedPadded(magnitude_of_symbolsPadded, 10,4)
function dopplerEstimationCentroidCorrectedPadded(magnitude_of_symbols, window_size, padding_factor)
    S = 0
    half_window = (window_size - 1) ÷ 2
    S_demodulated = argmax(magnitude_of_symbols)    # location of demodulated symbol

    # cuts out the surrounding of the magnitude to calculate centroid
    padded_magnitude_of_symbols, left_pad, right_pad = padded_window(magnitude_of_symbols, window_size, S_demodulated)
    values = padded_magnitude_of_symbols[S_demodulated-half_window+left_pad:S_demodulated+half_window+left_pad]
    n = length(values)    # Number of elements
  
    # Calculate the center of gravity
    center_of_gravity = sum(i * values[i] for i in 1:n) / sum(values) -1
    center_of_gravity = center_of_gravity 
    S_demodulated_centroid = (S_demodulated -1 - half_window + center_of_gravity) / padding_factor
    bin = Bw / (2^SF) # bin width in Hz
    ΔS1 = S_demodulated_centroid - S # symbol distance 
    ΔS2 = 2^SF - S_demodulated_centroid - S 
    if (ΔS1 < ΔS2)
        ΔS = ΔS1
    else
        ΔS = -ΔS2
    end
    fd_est = ΔS * bin # doppler estimation

    return S_demodulated_centroid, fd_est #, fd_error

end




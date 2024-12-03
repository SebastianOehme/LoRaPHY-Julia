function calcCentroid(magnitude_of_symbols)

    window_size = 10
    half_window = (window_size - 1) ÷ 2

    S_demodulated = argmax(magnitude_of_symbols)   # location of demodulated symbol

    # cuts out the surrounding of the magnitude to calculate centroid
    padded_magnitude_of_symbols, left_pad, right_pad = padded_window(magnitude_of_symbols, window_size, S_demodulated)
  
    values = padded_magnitude_of_symbols[S_demodulated-half_window+left_pad:S_demodulated+half_window+left_pad]
    n = length(values)    # Number of elements

    # Calculate the center of gravity
    center_of_gravity = sum(i * values[i] for i in 1:n) / sum(values) -1
    center_of_gravity = center_of_gravity 
    S_demodulated_centroid = S_demodulated - half_window + center_of_gravity
    
    return S_demodulated_centroid
end

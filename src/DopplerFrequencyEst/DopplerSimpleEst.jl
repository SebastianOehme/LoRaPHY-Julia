# simple arg max analysis
# TODO Error if Doppler is grater that bin * 2^SF 
# hard to tell if it wraped around already 


#magnitude_of_symbols_simple = performFFT(dechirp(chirp_signals[1]))
#plot(magnitude_of_symbols_simple)

#dopplerEstimationSimple(magnitude_of_symbols_simple)
#magnitude_of_symbols = magnitude_of_symbols_simple

function dopplerEstimationSimple(magnitude_of_symbols)
    S = 0
    S_demodulated = argmax(magnitude_of_symbols) - 1  # location of demodulated symbol
    bin = Bw / (2^SF) # bin width in Hz
    ΔS1 = S_demodulated - S 
    ΔS2 = 2^SF - S_demodulated - S 
    if (ΔS1 < ΔS2)
        ΔS = ΔS1
    else
        ΔS = -ΔS2
    end
    fd_est = ΔS  * bin # doppler estimation
    return fd_est
end




function performFFTPadded(signal::Vector{ComplexF64}, pad_length::Int)
    # Extend the signal with zeros to increase frequency resolution
    padded_signal = vcat(signal, zeros(ComplexF64, pad_length - length(signal)))
    
    fft_result = fft(padded_signal)
    magnitude_spectrum = abs.(fft_result)  # Magnitude of FFT

    return magnitude_spectrum
end
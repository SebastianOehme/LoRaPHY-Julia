using DSP

function dechirpWindowing(signal::Vector{ComplexF64})
    downchirp = ℯ.^(-1im * 2*pi*  (Bw/(2*T_s) .*t.^2 .+ -Bw/2 .*t))
    windowed_signal = signal .* hamming(length(signal)) # hanning()
    signalDechirped = windowed_signal .* downchirp

    return signalDechirped
end


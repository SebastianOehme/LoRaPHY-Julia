function dechirp(signal::Vector{ComplexF64})
    downchirp = ℯ.^(-1im * 2*pi*  (Bw/(2*T_s) .*t.^2 .+ -Bw/2 .*t))
    signalDechirped = signal .* downchirp

    return signalDechirped
end
using FFTW, Plots, Statistics, Noise


# ----------------- BASEBAND --------------------

# --------- Simplified Test ------------


# Parameters
#fs = 2*250000.0  # Sampling frequency (Hz)
#f0 = 125000.0     # Start frequency of the chirp (Hz)
#f1 = 250000.0    # End frequency of the chirp (Hz)
Δv = 8000   # [m/s]
c = 300000000 # lightspeed

Δv/c

# if LDRO is on SF = SF-2
SF_all = [7,8,9,10,11,12]
SF = SF_all[6]

fc_all = [433e6, 436e6, 868e6, 2100e6] # carrier in [Hz]
fc = fc_all[4]

Bw_all = [31.25, 62.5, 125, 250, 500].*1000
Bw = Bw_all[1]

S = 5
T_s =  (2^SF)/Bw # symbol periode
t = 0:1/(Bw):T_s  # Time vector from 0 to 0.1 seconds !!*100
f_offset = (Bw / (2^SF)) * S 
t_fold = (2^SF - S)/Bw 
stepsize = 1/(Bw) # !!* 100

# ------------- Generate chirp signal (frequency increases from f0 to f1 over time) --------------

t_to_tfold = 0:stepsize:t_fold
tfold_to_Ts =  t_fold+stepsize:stepsize:T_s

chirp_signal_1 = ℯ.^(1im * 2 * pi * ( Bw/(2*T_s) .* t_to_tfold.^2  .+ (f_offset - Bw/2) .* t_to_tfold ))
chirp_signal_2 = ℯ.^(1im * 2 * pi * ( Bw/(2*T_s) .* tfold_to_Ts.^2 .+ (f_offset - 3*Bw/2) .* tfold_to_Ts ))
chirp_signal = vcat(chirp_signal_1, chirp_signal_2)

#  ---------- Generate chirp signal with Doppler ------------

chirp_signal_doppler_1 = ℯ.^(1im * 2 * pi * ( (1 + Δv/c) .* (Bw/(2*T_s) .* t_to_tfold.^2 .+ (f_offset - Bw/2) * t_to_tfold) .+ (Δv/c .*fc .*t_to_tfold) ))
chirp_signal_doppler_2 = ℯ.^(1im * 2 * pi * ( (1 + Δv/c) .* (Bw/(2*T_s) .* tfold_to_Ts.^2 .+ (f_offset - 3*Bw/2) * tfold_to_Ts) .+ (Δv/c .*fc .*tfold_to_Ts) ))
chirp_signal_doppler = vcat(chirp_signal_doppler_1, chirp_signal_doppler_2)

#  ---------- Korrektor der Trägerfrequenz  --------
# mit den wissen über Δv und fc doppler der Trägerfrequenz ausrechnen
Δf = Δv/c*fc
correction = ℯ.^(- 1im * 2 * pi  .* Δf .* t)
chirp_signal_corrected = correction .* chirp_signal_doppler 

# ------------------ Test each step ------------------------

# ohne doppler ohne noise
plot(performFFT(dechirp(chirp_signal))[1:10])
# ohne doppler mit noise
plot!(performFFT(dechirp(add_noise(chirp_signal, -10)))[1:10])
# mit doppler ohne noise
plot!(performFFT(dechirp(chirp_signal_doppler))[1:45])
# mit doppler mit noise
plot!(performFFT(dechirp(add_noise(chirp_signal_doppler, 0)))[1:45])

# corrected ohne noise
plot!(performFFT(dechirp(chirp_signal_corrected))[1:10])
# corrected mit noise
plot!(performFFT(dechirp(add_noise(chirp_signal_corrected, 0)))[1:10])


# ------------------ Test each step ------------------------

# ------- Plot  Bandwidth / SER -----------

# ------- Plot  Δv / SER -----------

SER = Float64[]
SER_doppler = Float64[]

for Δv in 0:10:1000
    println(Δv)
    chirp_signal_doppler_1 = ℯ.^(1im * 2 * pi * ( (1 + Δv/c) .* (Bw/(2*T_s) .* t_to_tfold.^2 .+ (f_offset - Bw/2) * t_to_tfold) .+ (Δv/c .*fc .*t_to_tfold) ))
    chirp_signal_doppler_2 = ℯ.^(1im * 2 * pi * ( (1 + Δv/c) .* (Bw/(2*T_s) .* tfold_to_Ts.^2 .+ (f_offset - 3*Bw/2) * tfold_to_Ts) .+ (Δv/c .*fc .*tfold_to_Ts) ))
    chirp_signal_doppler = vcat(chirp_signal_doppler_1, chirp_signal_doppler_2)

    #SER = push!(SER, symbolErrorRate(chirp_signal, 0, 3000))
    SER_doppler = push!(SER_doppler, symbolErrorRate(chirp_signal_doppler, 0, 3000))

end

#plot(0:100:7500, SER, label="SER without Doppler", xlabel="SNR [dB]", ylabel="Symbol Error Rate", title="Comparison of SER of Non-/Dopplersignal", legend=:topright)
plot(0:10:1000, SER_doppler, label="SER with Doppler")


# ------- Plot SNR / SER -------------- 
SER = Float64[]
SER_doppler = Float64[]

for SNR_dB in -10:2:10
    println(SNR_dB)
    SER = push!(SER, symbolErrorRate(chirp_signal, SNR_dB, 100))
    SER_doppler = push!(SER_doppler, symbolErrorRate(chirp_signal_doppler, SNR_dB, 100))
end


plot(-10:2:10, SER, label="SER without Doppler", xlabel="SNR [dB]", ylabel="Symbol Error Rate", title="Comparison of SER of Non-/Dopplersignal", legend=:topright)
plot!(-10:2:10, SER_doppler, label="SER with Doppler")



function symbolErrorRate(signal::Vector{ComplexF64}, SNR_dB::Int64, loopPassages::Int64)
    
    symbol_cnt = 0

    for i in 1:1:loopPassages
        # add noise
        noisy_signal = add_noise(signal, SNR_dB)

        # Dechirping 
        sampled_signal = dechirp(noisy_signal)
        # FFT
        magnitude_spectrum = performFFT(sampled_signal)
        # Find Symbol
        symbol = argmax(magnitude_spectrum)-1

        if (symbol == S)
            symbol_cnt += 1
        end 
    end

    return symbol_error_rate = 1.0 - symbol_cnt/loopPassages
end




# ----------------- Dechirping -------------------

function dechirp(signal::Vector{ComplexF64})
    # generate downchirp (frequency decreases from f1 to f0 over time)
    # eigentlich nur ein minus ℯ.^(-1im ....) the complex conjugate of an upchirp is a downchirp
    downchirp = ℯ.^(-1im * 2*pi*  (Bw/(2*T_s) .*t.^2 .+ -Bw/2 .*t))
   
    # dechirp 
    signalDechirped = signal .* downchirp
    n = 1
    sampledSignal = signalDechirped[1:n:end-1]

    return sampledSignal
end


function performFFT(signal::Vector{ComplexF64})
    #N = length(signal)
    fft_result = fft(signal)
    #freqs = (0:N-1) * (fs / N)  # Frequency vector

    # Magnitude of FFT
    magnitude_spectrum = abs.(fft_result).^2

    # Plot FFT (magnitude vs frequency)
    #plot(freqs[1:div(N,2)], magnitude[1:div(N,2)], label="Magnitude", xlabel="Frequency (Hz)", ylabel="Magnitude", title="FFT of Chirp Signal", legend=:top)

    return magnitude_spectrum

end


# -------- Adding noise at a specific SNR -------- 

signal = fill(1.0 + 1.0im, 100)
SNR_dB = 0
plot!(1:1:100, sqrt.(real(signal).^2 + imag(signal).^2))
signal = add_noise(signal, 0)

signal = add_gauss(signal, 0.5)

function add_noise(signal::Vector{ComplexF64}, SNR_dB)

    # calculate the root mean square value
    signal_power = mean(abs2, signal)

    # from dB to linear
    SNR_linear = 10 ^ (SNR_dB/10)

    noise_power = signal_power / SNR_linear

    # generating complex gaussian noise
    noise = sqrt(noise_power / 2) * (randn(ComplexF64, length(signal)) + 1im * randn(ComplexF64, length(signal)))
    
    noisy_signal = signal + noise
    return noisy_signal

end

# ------------------- TEST NOISE FUNCTION -------------
signal = fill(1.0 + 1.0im, 100)

SNR_dB = -5
#function calculate_sigma_from_snr(signal, SNR_dB)
    # Step 1: Calculate signal power
    P_signal = mean(abs2, signal)
    
    # Step 2: Convert SNR from dB to linear scale
    SNR_linear = 10^(SNR_dB / 10)
    
    # Step 3: Calculate noise power
    P_noise = P_signal / SNR_linear
    
    # Step 4: Calculate the standard deviation of the noise
    sigma = sqrt(P_noise/2)
    
    # Step 5: gaussian
    signalNoise = add_gauss(signal, sigma)

    signalMyNoise = add_noise(signal, 0)
    plot(1:1:100, sqrt.(real(signalNoise).^2 + imag(signalNoise).^2))
    plot!(1:1:100, sqrt.(real(signal).^2 + imag(signal).^2))
    plot!(1:1:100, sqrt.(real(signalMyNoise).^2 + imag(signalMyNoise).^2))



# end

# Example usage
signal = fill(1 + 1im, 100)  # Example signal
SNR_dB = 20  # Desired SNR in dB

sigma = calculate_sigma_from_snr(signal, SNR_dB)
println("Calculated noise standard deviation: $sigma")





# ----- Plots -----
freq_n, magnitude_n = performFFT(noisy_signal)
plot(freq_n, magnitude_n, label="Noisy Chirp Signal", xlabel="Frequency (Hz)", ylabel="Magnitude", title="FFT of Chirp Signal", legend=:topright)

freq_n, magnitude_n =performFFT(dechirp(noisy_signal))
plot!(freq_n, magnitude_n, label="dechirped Doppler Signal", xlabel="Frequency (Hz)", ylabel="Magnitude", title="FFT of Chirp Signal", legend=:topright)











# -------- Calculation of the Bandwidth ----------

#freqs = freqs[1:div(N,2)]
#magnitude = magnitude[1:div(N,2)]

function bandwidthCalc(freqs::Vector{Float64}, magnitude::Vector{Float64})
    # Define a threshold to find significant frequencies
    threshold = 0.42 * maximum(magnitude)  # Set threshold as 50% of the peak

    # Find indices where the magnitude is greater than the threshold
    significant_indices = findall(magnitude .> threshold)

    # Find the frequencies corresponding to the significant indices
    significant_freqs = freqs[significant_indices]

    # Calculate bandwidth: difference between max and min significant frequencies
    bandwidth = maximum(significant_freqs) - minimum(significant_freqs)

    # Output the results
    println("Bandwidth of the signal: $bandwidth Hz")
    
end


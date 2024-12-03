# The idea is to have just one carrierFrequencyDopplerCorrection
# at the beginning of the packet which will result in errors


using FFTW, Plots, Statistics, Noise

# ------- LoRa specific parameters ------

# Spreading factor
SF_all = [7,8,9,10,11,12,18]
SF = SF_all[6]

# Carrier Frequency in [Hz]
fc_all = [433e6, 436e6, 868e6, 2100e6] 
fc = fc_all[1]

# Bandwidth in [Hz]
Bw_all = [31.25, 62.5, 125, 250, 500].*1000
Bw = Bw_all[5]

# Symbol to transmit 
S = [2,85,63,5,120,26,18]



# loop where every dependency on S has to be in
for q in 1:1:length(S)
    println(S[q])
end


# ------- LoRa specific calculated parameters ------
T_s =  (2^SF)/Bw             # Symbol periode in [s]
t = 0:1/(Bw):T_s -1/(Bw)           # Time vector from 0 to T_s
f_offset = (Bw / (2^SF)) * S[q] # starting frequency of the chirp in [Hz] 
t_fold = (2^SF - S)/Bw       # point in time when the frequency reaches the max and folds to min 
stepsize = 1/(Bw)            # sampling rate 


# ----  Constants -----
Δv = 1000     # relative velocity of receiver to transmitter [m/s]
c = 299792458 # speed of light in [m/s]

# ------------- Generate chirp signal (frequency increases from -Bw/2 to +Bw/2 over time) --------------
chirp_signal = Complex{Float64}[]

for t = 0:stepsize:T_s-1/(Bw)
    if (t <= t_fold)
        chirp_signal = push!(chirp_signal, ℯ^(1im * 2 * pi * ( Bw/(2*T_s) * t^2  + (f_offset - Bw/2) * t )))
    elseif (t > t_fold)
        chirp_signal = push!(chirp_signal,  ℯ^(1im * 2 * pi * ( Bw/(2*T_s) * t^2 + (f_offset - 3*Bw/2)* t )))
    end
end

#  ---------- Generate chirp signal with Doppler ------------
chirp_signal_doppler = Complex{Float64}[]

for t in 0:stepsize:T_s- 1/(Bw)
    if (t <= t_fold)
        chirp_signal_doppler = push!(chirp_signal_doppler,  ℯ^(1im * 2 * pi * ( (1 + Δv/c) * (Bw/(2*T_s) * t^2 + (f_offset - Bw/2) * t) + (Δv/c *fc *t) )))
    elseif (t > t_fold)
        chirp_signal_doppler = push!(chirp_signal_doppler, ℯ^(1im * 2 * pi * ( (1 + Δv/c) * (Bw/(2*T_s) * t^2 + (f_offset - 3*Bw/2) * t) + (Δv/c *fc *t) )))
    end
end


# ------------------ Add Noise ------------------
function add_noise(signal::Vector{ComplexF64}, SNR_dB)
    signal_power = mean(abs2, signal)       # calculate the root mean square value
    SNR_linear = 10 ^ (SNR_dB/10)           # from dB to linear
    noise_power = signal_power / SNR_linear # Calculate noise power
    noise = sqrt(noise_power / 2) * (randn(ComplexF64, length(signal)) + 1im * randn(ComplexF64, length(signal)))  # generating complex gaussian noise
    noisy_signal = signal + noise

    return noisy_signal
end


# ----------------- Dechirping ------------------

function dechirp(signal::Vector{ComplexF64})
    downchirp = ℯ.^(-1im * 2*pi*  (Bw/(2*T_s) .*t.^2 .+ -Bw/2 .*t))
    signalDechirped = signal .* downchirp

    return signalDechirped
end

# -------------------- FFT ----------------------
function performFFT(signal::Vector{ComplexF64})
    fft_result = fft(signal)
    magnitude_spectrum = abs.(fft_result)  # Magnitude of FFT

    return magnitude_spectrum
end

# --------------- Calculate SER -----------------
function symbolErrorRate(signal::Vector{ComplexF64}, SNR_dB::Int64, loopPassages::Int64)
    symbol_cnt = 0

    # to estimate the SER correctly the demodulation prozess has to be performed 
    # multiple (loopPassages-) times because the influence on the noise will be different every time 
    for i in 1:1:loopPassages
        noisy_signal = add_noise(signal, SNR_dB)         # add noise
        sampled_signal = dechirp(noisy_signal)           # Dechirping 
        magnitude_spectrum = performFFT(sampled_signal)  # FFT
        symbol = argmax(magnitude_spectrum) - 1          # Find Symbol

        # count the number of correct symbols
        if (symbol == S[q]) 
            symbol_cnt += 1
        end 
    end
    
    symbol_error_rate = 1.0 - symbol_cnt/loopPassages    # claculaton of SER
    return symbol_error_rate 
end


#  ---------- Correction of the carrier frequency  --------
function carrierFrequencyDopplerCorrection(signal::Vector{ComplexF64})
    Δf = Δv/c*fc                                 # calculating the frequency by which the carrier is doppler shifted
    correction = ℯ.^(- 1im * 2 * pi  .* Δf .* t) # complex signal with the frequency of -Δf 
    chirp_signal_corrected = correction .* signal # correct the doppler shift caused by the carrier 

    return chirp_signal_corrected
end


# --------------- Example Usecase  -----------------
# Comparison doppler signal vs non doppler signal vs doppler corrected signal

SER = Float64[]
SER_doppler = Float64[]
SER_doppler_corrected = Float64[]

for SNR_dB in -40:4:-20
    println(SNR_dB)
    SER = push!(SER, symbolErrorRate(chirp_signal, SNR_dB, 10000))
    #SER_doppler = push!(SER_doppler, symbolErrorRate(chirp_signal_doppler, SNR_dB, 10000))
    SER_doppler_corrected = push!(SER_doppler_corrected, symbolErrorRate(carrierFrequencyDopplerCorrection(chirp_signal_doppler), SNR_dB, 10000))
end


plot!(-40:4:-20, SER, label="signal without Doppler S=1", xlabel="SNR [dB]", ylabel="Symbol Error Rate", title="Comparison of SER of Non-/Doppler-/Corrected-Doppler-Signal \n SF= 7 Fc= 433MHz B= 500kHz Δv= 1000m/s",  titlefontsize = 10, legend=:bottomleft)
#plot!(-40:4:-20, SER_doppler, label="signal with Doppler")
plot!(-40:4:-20, SER_doppler_corrected, label="signal with corrected Dopplers S=1")
#savefig("C://Users//User//Desktop//Masterarbeit//DopplerPlots//SER_Comparison.png")



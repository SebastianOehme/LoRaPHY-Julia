
#function dopplerEstimationSimple(magnitude_of_symbols)
    
S_demodulated = argmax(magnitude_of_symbols) - 1  # location of demodulated symbol
bin = Bw / (2^SF) # bin width in Hz

# Peak Position Variant 1 
ΔS_positiv1 = S_demodulated  
ΔS_negativ1 = 2^SF - S_demodulated  

# Peak Position Variant 2
ΔS_positiv2 = S_demodulated  + 4096
ΔS_negativ2 = 2^SF - S_demodulated  + 4096

if (ΔS_positiv1 < ΔS_negativ1)
    ΔS = ΔS_positiv1
else
    ΔS = -ΔS_negativ1
end

# doppler estimation version 1 
fd_est_1 = ΔS_positiv2  * bin # doppler estimation

# doppler estimation version 2
fd_est_2 = ΔS  * bin


return fd_est
#end


t 

dopplerSignal_forCorrection_1 = ℯ.^(-1 .* im .* 2 .* π .*fd_est_1 .*  t )
dopplerSignal_forCorrection_2 = ℯ.^(-1 .* im .* 2 .* π .*fd_est_2 .*  t )

corrected_signal1 = chirp_signals[1] .* dopplerSignal_forCorrection_1
corrected_signal2 = chirp_signals[1] .* dopplerSignal_forCorrection_2

cor_signal_mag1 = performFFT(dechirp(corrected_signal1))
cor_signal_mag2 = performFFT(dechirp(corrected_signal2))

plot(cor_signal_mag1)
plot!(cor_signal_mag2)



#----------------------------------------------------
Bw = 31250.0 
T_s = 0.131072
nyquist_freq =  2*(50_000 + 100_000 + Bw/2 )
fs = 1_000_000 # Abtastrate in Hz um Nyquist zu bedienen fs​ > 2(Δf+f_IF​+2B​)

t_original = 0:1/Bw:T_s  # Zeitvektor
t_new = 0:1/fs:T_s
t_new = 0:1/fs:(T_s - 1/fs)   # Ensures correct number of samples


# Doppler korrektur factoren
fd_est_1 = 17662
fd_est_2 = 48912.048
fd_est_2 = -13587.952


# Empfangenes Signal im Basisband mit Doppler 4096-elemt Vector{Complex }
chirp_signals[1]

# Erzeugung eines signals für die Zwischenfrequenz
freq_erhoehung = 200_000
Signal_100KHz = ℯ.^(1 .* im .* 2 .* π .*freq_erhoehung .*  t_new )
mag  = performFFT(Signal_100KHz)
xfreqs = (fs / length(mag)) * (0:(length(mag) - 1))

plot(xfreqs[1:40000],mag[1:40000])

plot(real(Signal_100KHz[1:100]))

# chirp_signals muss geupsampled werden, sodass die Vectoren gleich lang sind 
upsample_factor = Int(round(fs / Bw))
chirp_signals_upsampled = resample(chirp_signals[1], upsample_factor)
plot(real(chirp_signals[1])[1:1000])
plot(real(chirp_signals_upsampled)[1:3200])

# Verschiebung auf Zwischenfrequenz
s_IF = chirp_signals_upsampled .* Signal_100KHz[1:length(chirp_signals_upsampled)]
pad_length = nextpow(2, length(s_IF))
padded_signal = vcat(s_IF, zeros(ComplexF64, pad_length - length(s_IF)))
fft_result1 = fft(padded_signal)
mag  = performFFT(padded_signal)
xfreqs = (fs / length(mag)) * (0:(length(mag) - 1))

plot(xfreqs[23000:31000],mag[23000:31000])

plot(real(s_IF)[1:100])

k = -Bw / T_s
downchirp = exp.(1im * π * k * t_new[1:length(chirp_signals_upsampled)] .^ 2)
downchirp = ℯ.^(-1im * 2*pi*  (Bw/(2*T_s) .*t_new[1:131056] .^2 .+ -Bw/2 .* t_new[1:131056]))


resample_factor = 1/32
resampled_signal = resample(padded_signal, resample_factor)
downchirp = ℯ.^(-1im * 2*pi*  (Bw/(2*T_s) .*t .^2 .+ -Bw/2 .* t))

signalDechirped1 = s_IF .* downchirp


padded_signal = vcat(signalDechirped1, zeros(ComplexF64, pad_length - length(signalDechirped1)))
fft_result1 = fft(signalDechirped1)
plot!(abs.(fft_result1))


mag  = performFFT(signalDechirped1)
plot(mag)
xfreqs = (fs / length(mag)) * (0:(length(mag) - 1))

plot(xfreqs[23000:31000],mag[23000:31000])

plot(xfreqs,mag)



resample_factor = 1/32
resampled_signal = resample(signalDechirped1, resample_factor)

pad_length = nextpow(2, length(resampled_signal))


padded_signal1 = vcat(resampled_signal, zeros(ComplexF64, pad_length - length(resampled_signal)))
fft_result1 = fft(padded_signal1)
magnitude_spectrum1 = abs.(fft_result1)
plot(magnitude_spectrum1)

padded_signal1 = vcat(signalDechirped1, zeros(ComplexF64, pad_length - length(signalDechirped1)))
fft_result1 = fft(padded_signal1)
magnitude_spectrum1 = abs.(fft_result1) 


xfreqs = (fs / length(fft_result1)) * (0:(length(fft_result1) - 1))
plot(xfreqs[1:40000], magnitude_spectrum1[1:40000])
argmax(magnitude_spectrum1)

# signal erzeugung für die frequenz korrektur
signalDechirped1
s_IF
dopplerSignal_forCorrection_1 = ℯ.^(-1 .* im .* 2 .* π .*fd_est_1 .*  t_new[1:131056] )
dopplerSignal_forCorrection_2 = ℯ.^(-1 .* im .* 2 .* π .*fd_est_2 .*  t_new[1:131056] )
plot(real(dopplerSignal_forCorrection_1)[1:100])

# Frequenz korrektur des verschoben signals
corrected_signal1 = s_IF .* dopplerSignal_forCorrection_1
corrected_signal2 = s_IF .* dopplerSignal_forCorrection_2

#corrected_signal1 = signalDechirped1 .* dopplerSignal_forCorrection_1

#dechirp 
downchirp = ℯ.^(-1im * 2*pi*  (Bw/(2*T_s) .*t_new[1:131056] .^2 .+ -Bw/2 .* t_new[1:131056]))

signalDechirped1 = corrected_signal1 .* downchirp
signalDechirped2 = corrected_signal2 .* downchirp

pad_length = 2^18
pad_length = nextpow(2, length(signalDechirped1))

padded_signal1 = vcat(signalDechirped1, zeros(ComplexF64, pad_length - length(signalDechirped1)))
padded_signal2 = vcat(signalDechirped2, zeros(ComplexF64, pad_length - length(signalDechirped2)))

fft_result1 = fft(padded_signal1)
magnitude_spectrum1 = abs.(fft_result1) 

fft_result2 = fft(padded_signal2)
magnitude_spectrum2 = abs.(fft_result2) 

# x axis
xfreqs = (fs / length(fft_result1)) * (0:(length(fft_result1) - 1))

plot(xfreqs[25000:40000], magnitude_spectrum1[25000:40000])
plot!(xfreqs[25000:40000], magnitude_spectrum2[25000:40000])
plot(xfreqs[20000:40000], magnitude_spectrum1[20000:40000])
plot!(xfreqs[18000:40000], magnitude_spectrum2[18000:40000])


1e5
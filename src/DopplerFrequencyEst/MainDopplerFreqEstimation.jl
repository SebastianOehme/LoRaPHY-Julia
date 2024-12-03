using FFTW, Plots, Statistics, Noise, LsqFit
include("../Utilities.jl")
include("DopplerCentroidEst.jl")
include("DopplerChirpCreation.jl")
include("DopplerSimpleEst.jl")
include("RecursiveLinearRegression.jl")
include("LinearRegression.jl")
include("PacketCreationRandom.jl")
include("ElevationAngle.jl")
#include("WholeOverpassFunction.jl")
include("CalcCentroid.jl")
include("LinearisationMedian.jl")
include("DechirpingWindowing.jl")



# ------- LoRa specific parameters ------

# Spreading factor
SF_all = [7,8,9,10,11,12, 16]
SF = SF_all[4]

# Carrier Frequency in [Hz]
fc_all = [433e6, 436e6, 868e6, 2100e6] 
fc = fc_all[1]

# Bandwidth in [Hz]
Bw_all = [31.25, 62.5, 125, 250, 500].*1000
Bw = Bw_all[1]

# Symbol to transmit 
number_of_packets = 1
packet_Length = 5158
preamble_length = 8

packet = fill_random_packet(packet_Length, preamble_length, SF)
packet = zeros(5128)

# ------- LoRa specific calculated parameters ------
T_s =  (2^SF)/Bw             # Symbol periode in [s]
t = 0:1/(Bw):T_s -1/(Bw)     # Time vector from 0 to T_s
stepsize = 1/(Bw)            # sampling rate 
fs_f = 1                     # oversampling factor

# ----  Constants -----
c = 299_792_458 # speed of light in [m/s]
H = 400_000; # orbit hight in meter
R = 6_371_000; # radius of earth
r = (H + R)/1000
g = 9.81;
G = 6.67430e-20  # Gravitational constant in m^3 kg^-1 s^-2
M_earth = 5.972e24  # Mass of the Earth in kg


# 1. Characterise the orbit and calculate the elevation angle

T_orbit  = 2 * pi * sqrt(r^3 / (G * M_earth))
t_velo = 0:1:2*T_orbit
elevationAngle = calcElevationAngle(t_velo)

# 2. Find the time periode of the first overpass in seconds
start_idx, end_idx = find_rise_fall_indices(elevationAngle)
t_overpass = start_idx:1:end_idx  
velo = calculationOfrealtiveVelocity(t_overpass)

plot(velo)
length(start_idx+510:T_s:end_idx)


# Comparison between constant velocity and non constant

chirp_signals_constant = Vector{Vector{Complex{Float64}}}(undef, length(start_idx:T_s:end_idx)) 
chirp_signals_nonconstant = Vector{Vector{Complex{Float64}}}(undef, length(start_idx:T_s:end_idx)) 

# chirp with constant Doppler


function constantDopplerchirp(start_idx, S)

    chirp_signal_doppler = Complex{Float64}[]

    f_offset = (Bw / (2^SF)) * S # starting frequency of the chirp in [Hz] 
    t_fold = (2^SF - S)/Bw       # point in time when the frequency reaches the max and folds to min 

    Δv = -calculationOfrealtiveVelocity(start_idx)

    for t in 0:stepsize:T_s- 1/(Bw)
        if (t <= t_fold)
            chirp_signal_doppler = push!(chirp_signal_doppler,  ℯ^(1im * 2 * pi * ( (1 + Δv/c) * (Bw/(2*T_s) * t^2 + (f_offset - Bw/2) * t) + (Δv/c *fc *t) )))
        elseif (t > t_fold)
            chirp_signal_doppler = push!(chirp_signal_doppler, ℯ^(1im * 2 * pi * ( (1 + Δv/c) * (Bw/(2*T_s) * t^2 + (f_offset - 3*Bw/2) * t) + (Δv/c *fc *t) )))
        end
    end
    
    return chirp_signal_doppler
end



for i in 0:1:length(start_idx:T_s:start_idx+packet_Length*T_s-T_s)-1
    println(i)
    chirp_signals_constant[Int(i)+1] = constantDopplerchirp(start_idx+i*T_s, packet[Int(i)+1])
    chirp_signals_nonconstant[Int(i)+1] = dopplerChirp(start_idx+i*T_s, packet[Int(i)+1])

end

# demodulation 

#magnitude_of_symbols_const = performFFTPadded(dechirpWindowing(chirp_signals_constant[1200]), 4*2^SF)
#magnitude_of_symbols_non_const = performFFTPadded(dechirpWindowing(chirp_signals_nonconstant[1200]), 4*2^SF)

magnitude_of_symbols_const = performFFT(dechirp(chirp_signals_constant[2200]))
magnitude_of_symbols_non_const = performFFT(dechirp(chirp_signals_nonconstant[2200]))

sampling_interval = 1  
time_in_seconds = (0:sampling_interval:(length(magnitude_of_symbols_const)-1)*sampling_interval) ./ 4

plot(time_in_seconds[7200:8000], magnitude_of_symbols_const[7200:8000])
plot!(time_in_seconds[7200:8000], magnitude_of_symbols_non_const[7200:8000])

plot(time_in_seconds[9000:9500], magnitude_of_symbols_const[9000:9500])
plot!(time_in_seconds[9000:9500], magnitude_of_symbols_non_const[9000:9500])


plot(magnitude_of_symbols_const)
plot!(magnitude_of_symbols_non_const)


plot(magnitude_of_symbols_const[300:350])
plot!(magnitude_of_symbols_non_const[300:350])

# ------------------------------- END ----------------------------


# 3.1 Calculating the chirp signal for a given Timestamp for whole Overpass

# Pre-allocate an array to store the results
chirp_signals = Vector{Vector{Complex{Float64}}}(undef, length(start_idx:T_s:end_idx))  # Assuming the loop runs 10 times


for i in 0:1:length(start_idx:T_s:start_idx+packet_Length*T_s-T_s)-1
    println(i)
    chirp_signals[Int(i)+1] = dopplerChirp(start_idx+i*T_s, packet[Int(i)+1])
end



# Add noise to the doppler distored chirp signals
noisyChirpSignals = Vector{Vector{Complex{Float64}}}(undef, length(start_idx:T_s:end_idx))  # Assuming the loop runs 10 times
SNR_dB = 0
for i in 0:1:length(start_idx:T_s:start_idx+packet_Length*T_s-T_s)-1
    println(i)
    noisyChirpSignals[Int(i)+1]  = add_noise(chirp_signals[Int(i)+1], SNR_dB)
end



# 3.2 Calculating the chirp signal for a given Timestamp just for a packet

# Pre-allocate an array to store the results at 
chirp_signals = Vector{Vector{Complex{Float64}}}(undef, length(start_idx+500:T_s:start_idx+500+T_s*packet_Length)-1)  


for i in 0:1:length(start_idx:T_s:start_idx+T_s*packet_Length-T_s)-1
    println(i)
    chirp_signals[Int(i)+1] = dopplerChirp(start_idx+500+i*T_s, packet[Int(i)+1])
end

# 4. Dechirp the signals 1 by 1 and perform an FFT

dopplerEstSimple = Vector{Float64}(undef, length(start_idx:T_s:start_idx+T_s*packet_Length)-1)
dopplerEstCentroid = Vector{Float64}(undef, length(start_idx:T_s:start_idx+T_s*packet_Length)-1)




# Estimate the Doppler from the preamble because here the Symbol values are known to the receiver

# for just one packet/ preamble
dopplerEstSimple = Vector{Float64}(undef, length(start_idx:T_s:start_idx+T_s*preamble_length)-1)
dopplerEstCentroid = Vector{Float64}(undef, length(start_idx:T_s:start_idx+T_s*preamble_length)-1)



length(start_idx:T_s:start_idx+T_s*preamble_length)
length(start_idx:T_s:start_idx+T_s*preamble_length)-1

for i in 1:1:length(start_idx:T_s:start_idx+T_s*preamble_length)-1
    println(i)
    #magnitude_of_symbols = performFFT(dechirp(chirp_signals[i]))
    magnitude_of_symbols = performFFTPadded(dechirpWindowing(chirp_signals[i]), 4*2^SF)
    #dopplerEstSimple[i] = dopplerEstimationSimple(magnitude_of_symbols)
    S_demodulated_centroid, dopplerEstCentroid[i] = dopplerEstimationCentroidCorrectedPadded(magnitude_of_symbols, 40, 4)
end


scatter!(dopplerEstCentroid)
plot(dopplerEstSimple)
plot(dopplerEstCentroid)



#------ loop -----
data_x = [1,2,3,4,5,6,7,8] # 1:8
data_x = 1:1:8
#data_x = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30] # 1:8
data_y = dopplerEstCentroid[1:8]

predictionIntervall = 1:1:preamble_length+packet_Length-8
predictionIntervall = 1:1:255
dopplerEstRecLinearMed, m, n  = liniearisationMedian(data_x, data_y, predictionIntervall)
dopplerEstRecLinearReg, m, n  = linearRegression(data_x, data_y, predictionIntervall)


for i in 1:255:5158-255
    data_y = dopplerEstCentroid[i:i+7]
    dopplerEstRecLinearReg, m, n  = linearRegression(data_x, data_y, predictionIntervall)
    plot!(i:254+i, dopplerEstRecLinearReg)
 
end


# Beispiel für den Frequenz Fehler eines Packets mit linear regression
data_y = dopplerEstCentroid[2000:2000+7]
dopplerEstRecLinearReg, m, n  = linearRegression(data_x, data_y, predictionIntervall)
plot(dopplerEstRecLinearReg)
dopplerEstRecLinearReg
plot!(dopplerEstCentroid[2000:2254])


scatter(dopplerEstRecLinearReg[1:10])
plot!(dopplerEstCentroid[2000:2010])

# Doppler estimation is done for the next symbols
S_true_est = zeros(Float64, packet_Length )
S_true =  zeros(Float64, packet_Length )
Symbol_error =  zeros(Float64, packet_Length )
Freq_error =  zeros(Float64, packet_Length )


data_y = dopplerEstCentroid[1:8]
predictionIntervall = 1:1:5128
dopplerEstRecLinearReg, m, n  = linearRegression(data_x, data_y, predictionIntervall)
plot(dopplerEstRecLinearReg[1:10])
plot!(dopplerEstCentroid[1:8])


m
n 
i = 181

# Newest version 15:53
# nochmal schauen warum sich der Geschätzte Doppler nicht ins negativ Biegt!
cdc = dopplerEstRecLinearReg
plot!(cdc[184:end])
plot(dopplerEstRecLinearReg_new)


last_x = 8
bin = Bw / (2^SF)



println(S_true)


plot(Symbol_error[2300:2350])
plot(diff(Symbol_error[1:1000]))

argmax(diff(Symbol_error[1:1000]))
plot(dopplerEstRecLinearReg[1:100])
plot!(dopplerEstCentroid[2400:2600])
minimumSum = sum(abs.((S_true - packet)[1:70]))
plot((S_true - packet)[1:100])
packet
S_true[1:10] - packet[1:10]







#for i in preamble_length+symbolSteps:symbolSteps:packet_Length


#dopplerEstRecLinear, m, n = linearRegression(data_x, data_y, predictionIntervall)
dopplerEstRecLinear, m, n  = liniearisationMedian(data_x, data_y, predictionIntervall)

# Doppler estimation is done for the next symbols
S_true_est = zeros(Float64, packet_Length )
S_true =  zeros(Float64, packet_Length )
Symbol_error =  zeros(Float64, packet_Length )
Freq_error =  zeros(Float64, packet_Length )







# Start Important Plot Centroid !!!
magnitude_of_symbols = performFFT(dechirp(chirp_signals[2790]))
plot(magnitude_of_symbols[1:20], label="data",  xlabel="symbols", ylabel="magnitude per symbol")
plot!(8:1:11, magnitude_of_symbols[12-4:11], label="window data points", seriestype = :scatter, mc=:red, ms=3, ma=1, xlabel="symbols", ylabel="magnitude per symbol", legendfontpointsize=7)
plot!(13:1:16, magnitude_of_symbols[13:12+4], label=false,seriestype = :scatter, mc=:red, ms=3, ma=1, xlabel="symbols", ylabel="magnitude per symbol")
plot!(12:1:12, magnitude_of_symbols[12:12], label="maximum",seriestype = :scatter, markershape=:cross, mc=:red, ms=4, ma=1, xlabel="symbols", ylabel="magnitude per symbol")


values = magnitude_of_symbols[12-4:12+4]
n = 9
center_of_gravity = sum(i * values[i] for i in 1:n) / sum(values) -1
x_v = center_of_gravity + 8
y_v = magnitude_of_symbols[12]

scatter!([x_v], [y_v], color = "green", label = "centroid", markersize = 3)
#savefig("C://Users//User//Desktop//Masterarbeit//DopplerPlots//CentroidCalc.png")
## End Important Plot !!!





# plot in seconds
x_values = (0:length(dopplerEstSimple)-1) .* T_s  # Generate x values in seconds
plot(x_values, dopplerEstSimple)


# Start Important Plot maximum vs centorid !!!
plot(dopplerEstSimple[1000:1100], label="maximum value estimation",  xlabel="time in [symbol periods]", ylabel="Doppler frequency in [Hz]")
plot!(dopplerEstCentroid[1000:1100], label="centroid estimation")
#savefig("C://Users//User//Desktop//Masterarbeit//DopplerPlots//maximum vs centorid.png")
# END Important Plot !!!



# ------------------ Add Noise ------------------
function add_noise(signal::Vector{ComplexF64}, SNR_dB)
    signal_power = mean(abs2, signal)       # calculate the root mean square value
    SNR_linear = 10 ^ (SNR_dB/10)           # from dB to linear
    noise_power = signal_power / SNR_linear # Calculate noise power
    noise = sqrt(noise_power / 2) * (randn(ComplexF64, length(signal)) + 1im * randn(ComplexF64, length(signal)))  # generating complex gaussian noise
    noisy_signal = signal + noise

    return noisy_signal
end

#randn(ComplexF64)
# checken 
noise_power = mean(abs2, noise)

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
        if (symbol == S) 
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

for SNR_dB in -44:1:-34
    println(SNR_dB)
    SER = push!(SER, symbolErrorRate(chirp_signal, SNR_dB, 1000))
    #SER_doppler = push!(SER_doppler, symbolErrorRate(chirp_signal_doppler, SNR_dB, 1000))
    SER_doppler_corrected = push!(SER_doppler_corrected, symbolErrorRate(carrierFrequencyDopplerCorrection(chirp_signal_doppler), SNR_dB, 1000))
end


plot(-44:1:-34, SER, label="signal without Doppler", xlabel="SNR [dB]", ylabel="Symbol Error Rate", title="Comparison of SER of Non-/Corrected-Doppler-Signal \n SF= 16 Fc= 433MHz B= 31.25kHz Δv= 7500m/s",  titlefontsize = 10, legend=:bottomleft)
#plot!(-40:4:-20, SER_doppler, label="signal with Doppler")
plot!(-44:1:-34, SER_doppler_corrected, label="signal with corrected Dopplers")
savefig("C://Users//User//Desktop//Masterarbeit//DopplerPlots//SER_Comparison_SF16.png")


println(SER)
println(SER - SER_doppler_corrected)
#x1 = SER - SER_doppler_corrected
plot(-44:1:-34, x1 , label="Durchgang 1", xlabel="SNR [dB]", ylabel="Difference of Symbol Error Rate", title="Difference of SER Between \n Non-/ and Corrected-Doppler-Signal") 
plot(-44:1:-34, SER - SER_doppler_corrected, label="Durchgang 2", xlabel="SNR [dB]", ylabel="Difference of Symbol Error Rate", title="Difference of SER Between \n Non-/ and Corrected-Doppler-Signal") 
savefig("C://Users//User//Desktop//Masterarbeit//DopplerPlots//SER_Comparison_Diff.png")




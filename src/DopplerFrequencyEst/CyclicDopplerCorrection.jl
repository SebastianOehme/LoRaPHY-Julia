using FFTW, Plots, Statistics, Noise, LsqFit, CircularArrayBuffers
include("../Utilities.jl")
include("DopplerCentroidEst.jl")
include("DopplerChirpCreation.jl")
include("DopplerGaussianEst.jl")
include("DopplerSimpleEst.jl")
include("LinearRegression.jl")
include("PacketCreationRandom.jl")
include("ElevationAngle.jl")
#include("WholeOverpassFunction.jl")
include("CalcCentroid.jl")
include("LinearisationMedian.jl")
include("DechirpingWindowing.jl")
include("FFTPadded.jl")



# ------- LoRa specific parameters ------

# Spreading factor
SF_all = [7,8,9,10,11,12, 16]
SF = SF_all[6]

# Carrier Frequency in [Hz]
fc_all = [433e6, 436e6, 868e6, 2100e6, 2400e6] 
fc = fc_all[5]

# Bandwidth in [Hz]
Bw_all = [31.25, 62.5, 125, 250, 500].*1000
Bw = Bw_all[1]


# ------- LoRa specific calculated parameters ------
T_s =  (2^SF)/Bw             # Symbol periode in [s]
t = 0:1/(Bw):T_s -1/(Bw)     # Time vector from 0 to T_s
stepsize = 1/(Bw)            # sampling rate 
bin = Bw / (2^SF)

# ----  Definition of the Satellite Orbit -----
c = 299_792_458; # speed of light in [m/s]
H = 400_000; # orbit hight in meter
R = 6_371_000; # radius of earth in meter
r = (H + R)/1000 # radius of satellite to earth center in km
g = 9.81;   # gravitational acceleration on earth
G = 6.67430e-20; # Gravitational constant in m^3 kg^-1 s^-2
M_earth = 5.972e24;  # Mass of the Earth in kg


# 1. Characterisation of the orbit and calculation of the elevation angle for two Orbits
T_orbit  = 2 * pi * sqrt(r^3 / (G * M_earth))
t_velo = 0:1:2*T_orbit
elevationAngle = calcElevationAngle(t_velo)

# 2. Find the time periode of the first overpass in seconds
start_idx, end_idx = find_rise_fall_indices(elevationAngle)
t_overpass = start_idx:1:end_idx  
#velo = calculationOfrealtiveVelocity(t_overpass)
#plot(velo)

end_idx - start_idx 
# Data to transmit 
T_s
number_of_packets = 20
packet_Length = 1000
preamble_length = 8

# for tests 
packet_Length = Int(floor((end_idx - start_idx)/T_s))
preamble_length = 8


# timestamp to start the transmission []

packet = fill_random_packet(packet_Length, preamble_length, SF)


# 3.1 Calculating the chirp signal for a given Timestamp for whole Overpass
# Pre-allocate an array to store the results
chirp_signals = Vector{Vector{Complex{Float64}}}(undef, length(start_idx:T_s:end_idx))  
# specific timestamp
start_idx = start_idx +200



for i in 0:1:length(start_idx:T_s:start_idx+packet_Length*T_s-T_s)-1
    println(i)
    chirp_signals[Int(i)+1] = dopplerChirp(start_idx+i*T_s, packet[Int(i)+1])
end


# Add noise to the doppler distored chirp signals
noisyChirpSignals = Vector{Vector{Complex{Float64}}}(undef, length(start_idx:T_s:end_idx))  
SNR_dB = -20
for i in 0:1:length(start_idx:T_s:start_idx+packet_Length*T_s-T_s)-1
    println(i)
    noisyChirpSignals[Int(i)+1]  = add_noise(chirp_signals[Int(i)+1], SNR_dB)
end



# 4. Dechirp the signals 1 by 1 and perform an FFT

dopplerEstSimple = Vector{Float64}(undef, length(start_idx:T_s:start_idx+T_s*packet_Length)-1)
dopplerEstCentroid = Vector{Float64}(undef, length(start_idx:T_s:start_idx+T_s*packet_Length)-1)


# Estimate the Doppler from the preamble because here the Symbol values are known to the receiver

# for just one packet/ preamble
dopplerEstSimple = Vector{Float64}(undef, length(start_idx:T_s:start_idx+T_s*preamble_length)-1)
dopplerEstCentroid = Vector{Float64}(undef, length(start_idx:T_s:start_idx+T_s*preamble_length)-1)



# Fine estimation for Doppler shift PLOT

magnitude1 = performFFT(dechirp(chirp_signals[1000]))
magnitude2 = performFFT(dechirpWindowing(chirp_signals[1000]))
magnitude3 = performFFTPadded(dechirpWindowing(chirp_signals[1000]), 4*2^SF)
plot(3270:3350, magnitude1[3270:3350], label="no window", xlabel="symbols", ylabel="magnitude per symbol")
plot(3270:3350,magnitude2[3270:3350], label=false)
scatter!(3270:3350,magnitude2[3270:3350], label="Without zero-padding", xlabel="symbols", ylabel="magnitude per symbol", color=1)
#savefig("C://Users//User//Desktop//Masterarbeit//MyMethodPlots//withoutpadding.png")

sampling_interval = 1  # Assuming each sample is 1 ms apart
# Create a time vector in seconds
time_in_seconds = (0:sampling_interval_ms:(length(magnitude3)-1)*sampling_interval_ms) ./ 4

scatter(time_in_seconds[3270*4:3350*4], magnitude3[3270*4:3350*4], label="with 4x zero-padding", xlabel="symbols", ylabel="magnitude per symbol")
#savefig("C://Users//User//Desktop//Masterarbeit//MyMethodPlots//withpadding.png")

#savefig("C://Users//User//Desktop//Masterarbeit//MyMethodPlots//hanning.png")






# Plots to explain the Novel Method 

# Big picture plot
magnitude_of_symbols_simple = performFFT(dechirp(chirp_signals[i]))

# Shift the plot by 40 x-values
shift_amount = -34
shifted_magnitude = circshift(magnitude_of_symbols_simple, shift_amount)

# Plot the original and shifted data
plot(magnitude_of_symbols_simple, label="S_Doppler", xlabel="symbols", ylabel="magnitude per symbol")
plot!(shifted_magnitude, label="S_true_estimate")
#savefig("C://Users//User//Desktop//Masterarbeit//MyMethodPlots//Bigpicture.png")

# high resolution plot
magnitude_of_symbols = performFFTPadded(dechirpWindowing(chirp_signals[1]), 4*2^SF)
# Shift the plot by 40 x-values
shift_amount = -34
shifted_magnitude = circshift(magnitude_of_symbols, shift_amount)

plot(magnitude_of_symbols[91:219],  label="S_Doppler", xlabel="symbols", ylabel="magnitude per symbol")
plot!(shifted_magnitude[91:219],  label="S_true_estimate")
#savefig("C://Users//User//Desktop//Masterarbeit//MyMethodPlots//Bigpicture2.png")

# Define the x-axis values, scaling to range 0 to 128
x_values = range(0, stop=128, length=length(magnitude_of_symbols))

argmax(magnitude_of_symbols)
magnitude_of_symbols[166]

# Plot with custom x-axis
plot(x_values[160:175], magnitude_of_symbols[160:175], label = "S_true_estimate", color = :2)
plot!([41.33], [magnitude_of_symbols[166]], seriestype = :scatter,  mc=:red, ms=3, ma=1, label = "max value")
# Create the array of length 128
data = zeros(128)
data[2] = 126  # Set the last value to 127
# Create the x values
x = 40:43
# Add a vertical line from the x-axis to the point at position 127
plot!([41, 41], [0, 126],xlabel="symbols", ylabel="magnitude per symbol", seriestype = :path, color = :black, lw = 1, label = "")
# Plot the points (discrete plot)
plot!(x, data, seriestype = :scatter,  mc=:blue, ms=3, ma=1, label = "S_true")
#savefig("C://Users//User//Desktop//Masterarbeit//MyMethodPlots//HighResolution.png")

#----------------------------------

for i in 1:1:length(start_idx:T_s:start_idx+T_s*preamble_length)-1
    println(i)
    #magnitude_of_symbols_simple = performFFT(dechirp(chirp_signals[i]))
    magnitude_of_symbols = performFFTPadded(dechirpWindowing(chirp_signals[i]), 4*4096)
    #dopplerEstSimple[i] = dopplerEstimationSimple(magnitude_of_symbols_simple)
    S_demodulated_centroid, dopplerEstCentroid[i] = dopplerEstimationCentroidCorrectedPadded(magnitude_of_symbols, 40, 4)
end

scatter(dopplerEstCentroid)
plot!(dopplerEstSimple)
plot(dopplerEstCentroid[1:end-2])
plot(diff(dopplerEstCentroid[1:end-2]))
plot(dopplerEstCentroid)

# workaround FFT wrapping 
x_cut = 0
difference = diff(dopplerEstCentroid)
y_value = 2*abs(dopplerEstCentroid[1])

for i in 1:1:length(dopplerEstCentroid)-2
    if abs(difference[i]) > Bw*0.9
        dopplerEstCentroid[1:i] .+= y_value
    end
end
plot(dopplerEstCentroid)


#------ loop -----

data_x = 1:1:8
data_y = dopplerEstCentroid[1:8]
circularbuffer = CircularArrayBuffer{Float64}(8)

for i in 1:1:length(data_y)
    push!(circularbuffer, data_y[i])
end


predictionIntervall = 1:1:9
dopplerEstRecLinearReg, m, n  = linearRegression(data_x, data_y, predictionIntervall)


scatter(data_x,data_y, label="preamble symbols")
plot!(dopplerEstRecLinearReg, label="1. linear regression",  xlabel="time in [T_s]", ylabel="Doppler shift [Hz]")
scatter!([9], [dopplerEstRecLinearReg[end]],  label="linear regression estimate for next symbol")
#savefig("C://Users//User//Desktop//Masterarbeit//MyMethodPlots//preambleEstimation.png")

push!(circularbuffer, 10084.674)
data_y = Vector{Float64}(circularbuffer)

dopplerEstRecLinearReg, m, n  = linearRegression(data_x, data_y, predictionIntervall)
scatter!([9], [10084.674],  label="next symbol")
plot!([2:10], dopplerEstRecLinearReg, label="linear regression update")
#savefig("C://Users//User//Desktop//Masterarbeit//MyMethodPlots//preambleEstimationNovel.png")




# Doppler estimation is done for the next symbols
S_true_est = zeros(Float64, packet_Length )
S_true =  zeros(Float64, packet_Length )
Symbol_error =  zeros(Float64, packet_Length )
Freq_error =  zeros(Float64, packet_Length )

i = 9
for i = preamble_length+1:1:packet_Length
    
    #i += 1
    # get a symbol outside the preamble
    magnitude_of_symbols = performFFTPadded(dechirpWindowing(chirp_signals[i]), 4*2^SF)

    # ------- PLOTS ------
    #plot(magnitude_of_symbols)

    # calculate its centroid 
    S_doppler, freq = dopplerEstimationCentroidCorrectedPadded(magnitude_of_symbols, 40, 4)
  
    # get symbols estimated doppler
    fd_est = dopplerEstRecLinearReg[end] # + bin

    # use this est doppler to calculate the real transmitted symbol
    ΔS = fd_est / bin

    # ------- PLOTS ------
    #plot(magnitude_of_symbols - ΔS)
    
    S_true_est[i] = ( S_doppler - ΔS ) % (2^SF) 
    if (S_true_est[i] < -1)
        S_true_est[i] = (2^SF) + S_true_est[i]
    end

    S_true[i] = round(Int, S_true_est[i])
    packet[i]

    # ------- PLOTS ------
    #plot(S_true_est)
    #plot!(S_true)

    Symbol_error[i] = S_true_est[i] -  S_true[i] 
    Freq_error[i] = Symbol_error[i] * bin

    freq = fd_est + Freq_error[i]
    #scatter!([10], [freq] )
  
    push!(circularbuffer, freq)
    data_y = Vector{Float64}(circularbuffer)

    dopplerEstRecLinearReg, m, n  = linearRegression(data_x, data_y, predictionIntervall)

end


plot(Freq_error)
plot((S_true - packet)[1:300])
println(S_true)
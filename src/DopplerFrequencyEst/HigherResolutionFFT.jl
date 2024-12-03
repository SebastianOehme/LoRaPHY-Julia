include("FFTPadded.jl")
include("DechirpingWindowing.jl")
include("DopplerChirpCreation.jl")
include("DopplerCentroidEst.jl")

packet_Length = 5158
preamble_length = 5158
padding_factor = 4

packet = fill_random_packet(packet_Length, preamble_length, SF)

t = 0:1/(Bw):T_s -1/(Bw)
stepsize = 1/(Bw)  

chirp_signals = Vector{Vector{Complex{Float64}}}(undef, length(start_idx:T_s:end_idx)) 

for i in 0:1:length(start_idx:T_s:start_idx+packet_Length*T_s-T_s)-1
    println(i)
    chirp_signals[Int(i)+1] = dopplerChirp(start_idx+i*T_s, packet[Int(i)+1])
end


dopplerEstSimple = Vector{Float64}(undef, length(start_idx:T_s:start_idx+T_s*packet_Length)-1)
dopplerEstCentroid = Vector{Float64}(undef, length(start_idx:T_s:start_idx+T_s*packet_Length)-1)
dopplerEstCentroidPaddedFFT = Vector{Float64}(undef, length(start_idx:T_s:start_idx+T_s*packet_Length)-1)

dopplerEstCentroid40 = Vector{Float64}(undef, length(start_idx:T_s:start_idx+T_s*packet_Length)-1)
dopplerEstCentroidPaddedFFT40 = Vector{Float64}(undef, length(start_idx:T_s:start_idx+T_s*packet_Length)-1)

i = 1
for i in 1:1:length(start_idx:T_s:start_idx+T_s*packet_Length-T_s)-1
    println(i)
    magnitude_of_symbolsPadded = performFFTPadded(dechirpWindowing(chirp_signals[i]),  2^SF*padding_factor)
    magnitude_of_symbols = performFFT(dechirpWindowing(chirp_signals[i]))


    dopplerEstSimple[i] = dopplerEstimationSimple(magnitude_of_symbols)
    #dopplerEstCentroid[i] = dopplerEstimationCentroidCorrectedPadded(magnitude_of_symbols, 10, 1)
    #dopplerEstCentroidPaddedFFT[i] = dopplerEstimationCentroidCorrectedPadded(magnitude_of_symbolsPadded, 10, padding_factor)

   # dopplerEstCentroid40[i] = dopplerEstimationCentroidCorrectedPadded(magnitude_of_symbols, 40, 1)
   S_demodulated_centroid, dopplerEstCentroidPaddedFFT40[i] = dopplerEstimationCentroidCorrectedPadded(magnitude_of_symbolsPadded, 40, padding_factor)
end


# ----- Comparison of real doppler Freq to estimated one with Packet full of S = 0 -----

t_overpass = start_idx+T_s:T_s:end_idx+T_s  
Δv = calculationOfrealtiveVelocity(t_overpass)

function calcDopplerShift(Δv::Vector{Float64})
    return f_doppler = (1 .+ -Δv./c) .* fc .- fc
end

f_doppler = calcDopplerShift(Δv)


# IMPORTANT FIGURE 
x_values = (0:length(f_doppler)-1) .* T_s
length(x_values)

plot(x_values[1:100], f_doppler[1:100], label="Caculated Doppler shift")
#plot!(x_values[300:500], dopplerEstSimple[300:500], label="Simple Doppler shift estimation", xlabel="time in [s]", ylabel="Doppler shift in [Hz]", )
plot(x_values[300:500], dopplerEstCentroid[300:500], label="Centroid Doppler shift estimation")
plot!(x_values[1:100], dopplerEstCentroidPaddedFFT40[1:100], label="Centroid Doppler shift estimation")

# high resolution Difference
plot(x_values[300:350], f_doppler[300:350], label="Caculated Doppler shift")
#plot!(x_values[300:500], dopplerEstSimple[300:500], label="Simple Doppler shift estimation", xlabel="time in [s]", ylabel="Doppler shift in [Hz]", )
plot!(x_values[300:350], dopplerEstCentroid[300:350], label="Centroid Doppler shift estimation")
plot!(x_values[300:350], dopplerEstCentroidPaddedFFT40[300:350], label="Centroid Doppler shift estimation")


#savefig("C://Users//User//Desktop//Masterarbeit//DopplerPlots//DopplerEstimationMulti.png")


error1 = f_doppler .- dopplerEstSimple
error2 = f_doppler .- dopplerEstCentroid
error7 = f_doppler .- dopplerEstCentroidPaddedFFT40

plot(x_values[1:5156], error1[1:5156], label="Simple Doppler shift estimation", xlabel="time in [s]", ylabel="error in [Hz]", legendfontpointsize=7)
plot!(x_values[1:5156], error2[1:5156], label="Centroid Doppler shift estimation\n window size = 10, FFT padding = 0")
plot!(x_values[1:5156], error7[1:5156], label="Centroid Doppler shift  estimation\n window size = 40, FFT padding = 4x")

#savefig("C://Users//User//Desktop//Masterarbeit//DopplerPlots//DopplerEstimationMultiError.png")



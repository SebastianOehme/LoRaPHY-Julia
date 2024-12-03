include("C://Users//User//Desktop//Masterarbeit//Simulation//src/LoRa_Demodulation/Dechirping.jl")
include("C://Users//User//Desktop//Masterarbeit//Simulation//Globals.jl")
using .Globals
using Plots
using ChirpSignal

# input should be the data symbols S = 0...2^SF-1
#=
data = [33]
test = modulationOfAPacket(data)

symbolsss = dechirpSignal(test)
println(symbolsss)
=#

function modulationOfAPacket(data::Int64, preamble_length::Int64, SF::Int64)

    preamble = Complex{Float64}[]
    frameSync = Complex{Float64}[]
    freqencySync = Complex{Float64}[]
    payload = Complex{Float64}[] 

    # upchirps for signal detection
    for i in 1:preamble_length
        preamble = vcat(preamble, chirpSignal(0, false))
    end


    # frame synchronc - to 0x0304 (1100000100) for public network in LoRaWAN 
    # TODO How to Encode the identifier above into the Chirps for know just random
    frameSync = vcat(frameSync, chirpSignal(rand(0:(2^SF)-1), false), chirpSignal(  rand(0:(2^SF)-1), false))


    # frequency sync - 2.25 base-downchirps
    #=
    for i in 1:3
        freqencySync = vcat(freqencySync, chirpSignal(0, true))
    end
    # adjust the length
    freqencySync225 = freqencySync[1:Int(round(length(freqencySync)  * 2.25/3))]
    =#

    # Chirping the data
    for j in 1:length(data)
        payload = vcat(payload, chirpSignal(data[j], false))
    end
   
    length(preamble)
    length(frameSync)
    #length(freqencySync225)
    length(payload)

    #LoRaPacket = vcat(preamble, frameSync, freqencySync225, payload)
    LoRaPacket = vcat(preamble, frameSync,  payload)

    
    x = 1:length(LoRaPacket)

    plot(x, real(LoRaPacket))

    return LoRaPacket


end


# -------------------------------

# Demodulation part for debugging 
#=

# sample signal with sample frequency
f_s = 1/params.BW

# create a base downchirp
DownChirp = chirpSignal(0, true)

S = Float64[]

Signal = LoRaPacket 


for i in length(DownChirp):length(DownChirp):length(Signal)
     
        # multiply the signal with a downchirp
        iszero = Signal[i-length(DownChirp)+1:i]  .* DownChirp;

        padded_vector = vcat(iszero, zeros(2^(params.SF+1) - 2^(params.SF)+1))
        # Perform the FFT
        fft_result = fft(padded_vector)

        # Calculate the magnitude of the FFT result
        magnitude_spectrum = abs.(fft_result)

     
        
        # calculate the symbol S
       
        S = push!(S, (argmax(magnitude_spectrum) - 1) /2)
       
end
println(S)

floorS = Int.(floor.(S))

println(floorS)

return S

=#

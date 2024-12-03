include("C://Users//User//Desktop//Masterarbeit//Simulation//Globals.jl")
include("C://Users//User//Desktop//Masterarbeit//Simulation//src/LoRa_Modulation/ChirpSignal.jl")
using .Globals
using Plots
using DSP
using FFTW


Signal = chirpSignal2_0(120, false)
println(dechirpSignal(Signal))

DownChirp = chirpSignal(0, true)
iszero = Signal  .* DownChirp;
plot(real(iszero))

function dechirpSignal(Signal::Vector{Complex{Float64}})

        # create a base downchirp
        DownChirp = chirpSignal(0, true)

        S = Float64[]
        
        for i in length(DownChirp):length(DownChirp):length(Signal)
             
                # multiply the signal with a downchirp
                iszero = Signal[i-length(DownChirp)+1:i]  .* DownChirp;

                # downsampling signal with BW
                # define the sampling frequency depends on the resolution/stepsize in ChirpSignal.jl
                # it should be 1/params.BW. In ChirpSignal.jl the stepsize is 1/(params.BW*100) therefore
                # n must be 100. 
                # additionally the vector is of length 2^SF + 1 and must be reduced to 2^SF
                # otherwise the FFT can not be performed in a sufficient way. Needs a vector of 2^x 
                # instead of cutting padding could be another solution 
                # padded_vector = vcat(iszero, zeros(2^(params.SF+1) - 2^(params.SF)+1))
                n = 100
                sampled_signal = iszero[1:n:end-1]
                
                # Perform the FFT
                fft_result = fft(sampled_signal)

                # Calculate the magnitude of the FFT result
                magnitude_spectrum = abs.(fft_result)

                # calculate the symbol S
                S = push!(S, argmax(magnitude_spectrum) -1)
        end

        return S
end


# -------------------- TESTS ---------------------



       


 

   



  
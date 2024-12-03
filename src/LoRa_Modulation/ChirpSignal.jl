include("C://Users//User//Desktop//Masterarbeit//Simulation//Globals.jl")
using .Globals
using Plots
using DSP


function chirpSignal(S, isdownchirp::Bool)


    Chirp = Complex{Float64}[]
    BaseChirp = Complex{Float64}[]

    f_offset = (params.BW / (2^params.SF)) * S 
    t_fold = (2^params.SF - S)/params.BW 
    T_s =  2^params.SF/params.BW # symbol periode
    # Simulations stepsize - the "sampling frequency" should be equal to bandwidth (chinesisches Paper)
    stepsize = 1/(params.BW*100)
    

    if (isdownchirp == 1)
        for t = 0:stepsize:T_s
            BaseChirp = push!(BaseChirp, ℯ^(1im * 2*pi*  (params.BW*t^2 / (2*T_s) + (-params.BW/2)*t)) )
        end
        # the complex conjugate of an upchirp is a downchirp
        Chirp = real(BaseChirp) - 1im * imag(BaseChirp)
    else
        for t = 0:stepsize:T_s
            if (t <= t_fold)
                Chirp = push!(Chirp, ℯ^(1im * 2*pi*  (params.BW*t^2 / (2*T_s) + (f_offset - params.BW/2)*t)) )
            elseif (t > t_fold)
                Chirp = push!(Chirp, ℯ^(1im * 2*pi*  (params.BW*t^2 / (2*T_s) + (f_offset - 3*params.BW/2)*t)))
            end
        end
    end


    return Chirp
    
end


function chirpSignal2_0(S, isdownchirp::Bool)

    Chirp = Complex{Float64}[]
    BaseChirp = Complex{Float64}[]

    f_offset = (params.BW / (2^params.SF)) * S 
    t_fold = (2^params.SF - S)/params.BW 
    T_s =  2^params.SF/params.BW # symbol periode
    # Simulations stepsize - the "sampling frequency" should be equal to bandwidth (chinesisches Paper)
    stepsize = 1/(params.BW*100)
    

    if (isdownchirp == 1)
        for t = 0:stepsize:T_s
            BaseChirp = push!(BaseChirp, ℯ^(-1im * 2*pi*  (params.BW*t^2 / (2*T_s) + (-params.BW/2)*t)) )
        end
    else
        for t = 0:stepsize:T_s
            if (t <= t_fold)
                Chirp = push!(Chirp, ℯ^(1im * 2*pi*  (params.BW*t^2 / (2*T_s) + (f_offset - params.BW/2)*t)) )
            elseif (t > t_fold)
                Chirp = push!(Chirp, ℯ^(1im * 2*pi*  (params.BW*t^2 / (2*T_s) + (f_offset - 3*params.BW/2)*t)))
            end
        end
    end


    return Chirp
    
end



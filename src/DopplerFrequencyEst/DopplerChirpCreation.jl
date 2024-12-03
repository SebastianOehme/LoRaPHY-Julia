#  ---------- Generate chirp signal with simulated Doppler Δv is calculated ------------

function dopplerChirp(start_idx, S)

    chirp_signal_doppler = Complex{Float64}[]

    f_offset = (Bw / (2^SF)) * S # starting frequency of the chirp in [Hz] 
    t_fold = (2^SF - S)/Bw       # point in time when the frequency reaches the max and folds to min 

    for t in 0:stepsize:T_s - 1/(fs_f*Bw) # achtung hier eigentlich - 1/(Bw)
        if (t <= t_fold)
            Δv = -calculationOfrealtiveVelocity(start_idx+t);

            chirp_signal_doppler = push!(chirp_signal_doppler,  ℯ^(1im * 2 * pi * ( (1 + Δv/c) * (Bw/(2*T_s) * t^2 + (f_offset - Bw/2) * t) + (Δv/c *fc *t) )));
        elseif (t > t_fold)
            Δv = -calculationOfrealtiveVelocity(start_idx+t);

            chirp_signal_doppler = push!(chirp_signal_doppler, ℯ^(1im * 2 * pi * ( (1 + Δv/c) * (Bw/(2*T_s) * t^2 + (f_offset - 3*Bw/2) * t) + (Δv/c *fc *t) )));
        end
    end

    return chirp_signal_doppler
end



#  ---------- relative velocity of a LEO satellite ------------

function calculationOfrealtiveVelocity(t_velo)
    # satellites angular position
    ϕ = t_velo .* sqrt.(g ./R) .* (1 .+ H ./R).^-1.5
    Δv = sqrt.(g.*R./(1 .+H./R)) .* sin.(ϕ) ./ sqrt.( (1 .+ H ./ R).^2  .- 2 .*(1 .+ H ./ R).*cos.(ϕ)  .+ 1)
    return Δv
end


function calculationOfrealtiveVelocityDiff(t_velo)
    ω = sqrt(g*R^2)  / (R+H)^1.5
    ω_t = t_velo .* ω 
    r_s = R + H
    Δv = ω .* r_s .* R .* sin.(ω_t) ./ (sqrt.(r_s^2 .+ R^2 .- 2 .* r_s .* R .* cos.(ω_t)))

    return Δv
end


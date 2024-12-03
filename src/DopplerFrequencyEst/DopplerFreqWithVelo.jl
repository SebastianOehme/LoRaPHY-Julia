fc
c

t_overpass = start_idx:1:end_idx  
Δv = calculationOfrealtiveVelocity(t_overpass)

plot(-Δv)

function calcDopplerShift(Δv::Vector{Float64})
    return f_doppler = (1 .+ -Δv./c) .* fc .- fc
end


plot(f_doppler)
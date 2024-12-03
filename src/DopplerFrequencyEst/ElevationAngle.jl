function calcElevationAngle(t_velo)
    β = 1 + H/R

    ω = sqrt(g / R) * β^(-1.5)  # Angular velocity
    ϕ = ω .* t_velo # Satellite angular position

    # Distance between satellite and ground station
    d = R .* sqrt.(β^2 .+ 1 .- 2 .* β .* cos.(ϕ))

    # Elevation angle in radians (with clamping to prevent DomainError Floating point)
    ε = asin.(clamp.((R .* (β .* cos.(ϕ) .- 1)) ./ d, -1.0, 1.0))

    # Convert elevation angle to degrees if needed
    ε_deg = ε .* (180 / π)
    return ε_deg 
end
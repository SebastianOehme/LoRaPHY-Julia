
##################################################################

# https://juliaspace.github.io/SatelliteToolbox.jl/stable/tutorials/iss_observation/

# earth rotation? 

##################################################################

# Zielstellung ist:
# 1. TLE auslesen und damit die Geschwindigkeit herausfinden
# 2. den Winkel berechnen
# 3. den Doppler shift berechnen 
# 4. den Doppler shift plotten... 

using SatelliteToolbox
using Dates
using PhysicalConstants.CODATA2018
using Unitful
using LinearAlgebra
using Plots

# Location of Würzburg
latitude_wue = 49.48; # [°]
longitude_wue = 9.56; # [°]
altitude_wue = 177; # [m] 

# Location of Frankfurt
latitude_fra = 50.7; # [°]
longitude_fra = 8.41; # [°]
altitude_fra = 177; # [m] 


F_C = 2100e6; # [Hz]
c_0 = 299792458; # [m/s^-1]
R_e = 6371000;
B = 31_250


# get latest TLE from i.e. Sonate
fetcher = create_tle_fetcher(CelestrakTleFetcher)
tles = fetch_tles(fetcher; satellite_name = "ISS (ZARYA)")

iss_tle = tles[1]

orbp = Propagators.init(Val(:SGP4), iss_tle) 


# propagation for 5400s or 90 min ~one orbit 
ret = Propagators.propagate!.(orbp, 0:1:86400)

# extract the position vector in earth center inertial (ECI)
position_teme = first.(ret)

#extract the velocity vector in earth center inertial (ECI)
velocity_teme = [subtuple[2] for subtuple in ret]


# vector r and velocity v pseudo earth fixed (p ECEF) The code Propagators.epoch(orbp) .+ (collect(0:1:86400) ./ 86400) obtains the Julian Day [UTC] of each propagation instant.
vr_pef = r_eci_to_ecef.(TEME(), PEF(), Propagators.epoch(orbp) .+ (collect(0:1:86400) ./ 86400)) .* position_teme
vv_pef = r_eci_to_ecef.(TEME(), PEF(), Propagators.epoch(orbp) .+ (collect(0:1:86400) ./ 86400)) .* velocity_teme



# -------- ploting  position in ECEF coordinate system --------
x_coords = [point[1] for point in vr_pef]
y_coords = [point[2] for point in vr_pef]
z_coords = [point[3] for point in vr_pef]
scatter3d(x_coords, y_coords, z_coords, xlabel="X Axis", ylabel="Y Axis", zlabel="Z Axis", title="3D Scatter Plot of Coordinates", markersize=1, legend=false)



velocity = [sqrt(sum(vec .^2)) for vec in vv_pef]

# converting IoT sensor position to North-East-Down coordinate system
vr_ned = ecef_to_ned.(vr_pef, latitude_wue |> deg2rad, longitude_wue |> deg2rad, altitude_wue; translate = true)


#compute distance from IoT sensor to satellite
distance_to_GS = [sqrt(sum(vec .^2)) for vec in vr_ned]
plot(0:1:85999,distance_to_GS[1:86000]/1000)

# idee 
# über die distance des Satelliten zur GS kann man
# auch seine relative Geschwindigkeit zur GS ausrechen
vel_to_GS = diff(distance_to_GS)
#vel_to_GS=abs.(vel_to_GS)
t = 1:1:86400
plot(0:1:9999, vel_to_GS[1:10000], label=[] )
xlabel!("t in [s]")
ylabel!("relative velocity in [m/s]")
#savefig("C://Users//User//Desktop//Masterarbeit//DopplerPlots//Velocity.png")



# jetzt die Geschwindigkeitsänderung pro sekunde zur GS
a_sat = diff(vel_to_GS)

t = 1:1:86399
plot(t, a_sat)
#plot(0:1:325, a_sat[77800:78125], label=[vλ])
title!(" ")
xlabel!("t in [s]")
ylabel!("acceleration in [m/s^2]")
#savefig("C://Users//User//Desktop//Masterarbeit//DopplerPlots//Acceleration.png")


# jetzt noch den Doppler ausrechnen
# -------------- rate
f_d1 =  (a_sat./300000000) .* 436.7e6
f_d2 =  (a_sat./300000000) .* 868e6
f_d3 =  (a_sat./300000000) .* 2100e6

plot(0:1:99, f_d1[1:100],   label="436.7 MHz")
plot!(0:1:99, f_d2[1:100], label="868 MHz")
plot!(0:1:99, f_d3[1:100], label="2.1 GHz")
xlabel!("t in [s]")
ylabel!("Doppler rate [Hz/s]")
#savefig("C://Users//User//Desktop//Masterarbeit//DopplerPlots//Dopplerrate.png")

# -------------- shift

f_d1 =  (((1 ./  (1 .+ (vel_to_GS./300000000))) .* 436.7e6) .- 436.7e6)./1000
f_d2 =  (((1 ./  (1 .+ (vel_to_GS./300000000))) .* 868e6) .- 868e6)./1000
f_d3 =  (((1 ./  (1 .+ (vel_to_GS./300000000))) .* 2100e6) .- 2100e6)./1000

plot(5000:1:7000, f_d1[5000:7000],   label="436.7 MHz")


plot!(0:1:9999, f_d2[1:10000], label="868 MHz")
plot!(0:1:9999, f_d3[1:10000], label="2.1 GHz")
xlabel!("t in [s]")
ylabel!("Doppler shift in [kHz]")
#savefig("C://Users//User//Desktop//Masterarbeit//DopplerPlots//Dopplershift.png")

f_dd1 = diff(f_d1)*1000
f_dd2 = diff(f_d2)*1000
f_dd3 = diff(f_d3)*1000
plot(0:1:325, f_dd1[77800:78125],   label="436.7 MHz")
plot!(0:1:325, f_dd2[77800:78125], label="868 MHz")
plot!(0:1:325, f_dd3[77800:78125], label="2.1 GHz")
xlabel!("t in [s]")
ylabel!("Doppler rate [Hz/s]")
#savefig("C://Users//User//Desktop//Masterarbeit//DopplerPlots//Dopplerrate.png")
# calculate angle beta (the angle between the groundstationToSatellite vector in NED 
# and the velocity vecotr of the satellite in ECEF)




vec_GS_to_Sat = vr_ned
vec_sat_velo = vv_pef;

# -------- ploting  velocity vector in ECEF --------
x_coords = [point[1] for point in vec_sat_velo]
y_coords = [point[2] for point in vec_sat_velo]
z_coords = [point[3] for point in vec_sat_velo]
scatter3d(x_coords, y_coords, z_coords, xlabel="X Axis", ylabel="Y Axis", zlabel="Z Axis", title="3D Scatter Plot of Coordinates", markersize=1, legend=false)





# calculating the velocity vector with delta position
#
#
# .......


#calculating angle between the two vectors
beta = acosd.(clamp.(vec_GS_to_Sat.⋅vec_sat_velo./(norm.(vec_GS_to_Sat).*norm.(vec_sat_velo)), -1, 1))

#plot(0:60:5400, beta)


# computation of the elevation
vλ = map(v -> atand(-v[3], sqrt(v[1]^2 + v[2]^2)), vr_ned)
plot(0:1:9999, vλ[1:10000], label=[vλ] )
#title!("")
xlabel!("t in [s]")
ylabel!("Elevation angle in [°]")
#savefig("C://Users//User//Desktop//Masterarbeit//DopplerPlots//Elevation.png")
# find elevation greater than 0
ids = findall(>=(0), vλ)
println(ids)
vλ[1:377]

function find_sequences(arr)
    # Initialisieren der Variablen
    sequences = []
    start_idx = 1
    
    while start_idx <= length(arr)
        end_idx = start_idx
        while end_idx < length(arr) && arr[end_idx + 1] == arr[end_idx] + 1
            end_idx += 1
        end
        
        if end_idx > start_idx
            push!(sequences, (start_idx, end_idx))
        end
        
        start_idx = end_idx + 1
    end
    
    return sequences
end


# Find the start and end of the overpasses
overpasses_start_end = find_sequences(ids)


overpasses_start_end[1][2]


#number of overpasses in the time intervall
length(overpasses_start_end)

# Ergebnis anzeigen
println(overpasses_start_end)

# 
#overpasses_start_end[2][2]


# -------------- sequence + steps ---------
function find_sequences_elements_and_steps(arr)
    # Initialize variables
    sequences = []
    total_steps = 0
    start_idx = 1
    
    while start_idx <= length(arr)
        end_idx = start_idx
        while end_idx < length(arr) && arr[end_idx + 1] == arr[end_idx] + 1
            end_idx += 1
        end
        
        if end_idx > start_idx
            push!(sequences, arr[start_idx:end_idx])
            # Count the steps in this sequence
            total_steps += (end_idx - start_idx)
        end
        
        start_idx = end_idx + 1
    end
    
    return sequences, total_steps
end



# Find sequences and count steps
sequences, total_steps = find_sequences_elements_and_steps(ids)

F_R = Float64[]  # Initialize an empty array for Hz

sequences[2][overpasses_start_end[2][2]]

#overpasses_start_end[2][2] - overpasses_start_end[2][1]

for i in 1:length(overpasses_start_end)
    println(i)
    end_of_overp = overpasses_start_end[i][2] - overpasses_start_end[i][1]
    println(end_of_overp)
    F_R_calc = 1 ./(1 .+ (velocity[sequences[i][1]:sequences[i][end_of_overp]]./c_0.* cos.(beta[sequences[i][1]:sequences[i][end_of_overp]].*pi./180))).* F_C
    F_R = append!(F_R, F_R_calc)
end



# received carrier frequency
#F_R = 1 ./(1 .+ (velocity[sequences[1][1]:sequences[1][106]]./c_0.* cos.(beta[sequences[1][1]:sequences[1][106]].*pi./180))).* F_C
F_R
F_D = F_R .- F_C

difFD = diff(F_D)
length(diff(F_D))

for i in 1:length(diff(F_D))
    if (difFD[i] > 5000)
        println(i+1)
    end
end

plot(F_D[377:1027-1])





plot(0:60:(total_steps-1)*60, F_D, linewidth=2)
title!("Doppler Shift from all Overpasses from 31.07.2024 ")
xlabel!("t in seconds")
ylabel!("Hz")

sequences[2][1]
#---------------------------------------------------------------------------------------------------------------

#find the exact dates of overpass
#times_of_overpass = getindex(Propagators.epoch(orbp) .+ (collect(0:60:86400) ./ 86400) .|> julian2datetime, ids)

#using write method
#CSV.write("data/times_of_overpass.csv", DataFrame(times_of_overpass), header = false)FF
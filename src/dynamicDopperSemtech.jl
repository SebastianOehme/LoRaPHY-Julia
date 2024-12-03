using Plots

# Spreading factor
SF_all = [7.0,8.0,9.0,10.0,11.0,12.0]
SF = SF_all[6]

# Carrier Frequency in [Hz]
fc_all = [433e6, 436e6, 868e6, 2100e6] 
fc = fc_all[1]

# Bandwidth in [Hz]
Bw_all = [31.25, 62.5, 125, 250, 500].*1000
Bw = Bw_all[1]


# ------- LoRa specific calculated parameters ------

F_dynamic =  Array{Float64}(undef, 5, 6)
T_sM =  Array{Float64}(undef, 5, 6)


for i in 1:1:5
    for j in 1:1:6
        T_s =  (2^SF_all[j])/Bw_all[i]
        T_sM[i,j] = T_s 
        if (T_s > 0.016)
            L = 16.0
        else
            L = 1.0
        end
        F_dynamic[i,j] = L * Bw_all[i] / (3 * 2^SF_all[j])
    end
end

println(round.(F_dynamic,digits=1))

T_sM
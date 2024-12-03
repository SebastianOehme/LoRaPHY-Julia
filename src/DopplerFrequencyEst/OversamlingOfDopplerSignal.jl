include("DopplerChirpCreation.jl")
include("DopplerSimpleEst.jl")
include("DechirpingWindowing.jl")
include("FFTPadded.jl")
using DSP

# ------- LoRa specific parameters ------

# Spreading factor
SF_all = [7,8,9,10,11,12, 16]
SF = SF_all[4]

# Carrier Frequency in [Hz]
fc_all = [433e6, 436e6, 868e6, 2100e6] 
fc = fc_all[4]

# Bandwidth in [Hz]
Bw_all = [31.25, 62.5, 125, 250, 500].*1000
Bw = Bw_all[1]





# sampling factor
fs_f = 2
# ------- LoRa specific calculated parameters ------
T_s = (2^SF)/Bw             # Symbol periode in [s]
t = 0:1/(fs_f*Bw):T_s -1/(fs_f*Bw)  # Time vector from 0 to T_s 
length(t)
stepsize = 1/(Bw)            # sampling rate 
bin = Bw / (2^SF)


stepsize = 1/(fs_f*Bw)          # over sampling rate

# ----  Constants -----
c = 299_792_458 # speed of light in [m/s]
H = 480_000; # orbit hight in meter
R = 6_371_000; # radius of earth
r = (H + R)/1000
g = 9.81;
G = 6.67430e-20  # Gravitational constant in m^3 kg^-1 s^-2
M_earth = 5.972e24  # Mass of the Earth in kg


# 1. Characterise the orbit and calculate the elevation angle

T_orbit  = 2 * pi * sqrt(r^3 / (G * M_earth))
t_velo = 0:1:2*T_orbit
elevationAngle = calcElevationAngle(t_velo)

# 2. Find the time periode of the first overpass in seconds
start_idx, end_idx = find_rise_fall_indices(elevationAngle)
t_overpass = start_idx:1:end_idx  

5309 - 5985
676/2
S = 100

# Get the ECB Doppler shifted signal
dopplerChirpSignal = dopplerChirp(5985-338, S)


# simplest approach 
dopplerDechirpSignal = dechirp(dopplerChirpSignal)
fft_result = fft(dopplerDechirpSignal)
plot(imag(fft_result))
magnitude_spectrum = abs.(fft_result) 
plot(magnitude_spectrum)

# -------------------- FFT ----------------------
 # get a symbol outside the preamble
 magnitude_of_symbols = performFFTPadded(dechirpWindowing(dopplerChirpSignal), 4*fs_f*2^SF)
 # calculate its centroid 
S_doppler, freq = dopplerEstimationCentroidCorrectedPadded(magnitude_of_symbols, 40, 4)


plot(magnitude_of_symbols)

frequencyOfDopplerSignal = performFFT(dopplerDechirpSignal)
plot(frequencyOfDopplerSignal)

#frequencyOfDopplerSignal = performFFT(dopplerDechirpSignal)
#S_doppler = argmax(frequencyOfDopplerSignal) - 1


# combine two peaks
# find the two maxima (idx and absvalue)

# Find indices of all local maxima
function find_local_maxima(data)
    return [i for i in 2:length(data)-1 if data[i-1] < data[i] && data[i] > data[i+1]]
end

# Find indices of the two largest peaks
function find_two_peaks(data)
    local_maxima = find_local_maxima(data)
    # Get magnitudes and sort the maxima by descending magnitude
    peak_magnitudes = data[local_maxima]
    sorted_indices = sortperm(peak_magnitudes, rev=true)
    # Get the top two peaks
    top_peaks = local_maxima[sorted_indices[1:2]]
    return top_peaks, data[top_peaks]
end

# Find the indices and magnitudes of the two largest peaks
peak_indices, peak_magnitudes = find_two_peaks(magnitude_of_symbols)

# for positve Doppler use the peak that is further away 
peak_indices[1]  
peak_indices[2] 

# for negative Doppler the peak that is closer

# use peak 1 to calculate Doppler
Δf1 = (( peak_indices[1] / 4 ) - S ) * bin

# use peak 2 to calculate Doppler
Δf2 = 16384/4* bin - (( peak_indices[2] / 4 ) - S ) * bin


Δf1 - Δf2




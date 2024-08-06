#package
import numpy as np
import pandas
from decimal import Decimal
from scipy.interpolate import UnivariateSpline
from scipy.signal import find_peaks

#scientific figure format
def sci_notation(number, sig_fig=2):
    ret_string = "{0:.{1:d}e}".format(number, sig_fig)
    siga, sigb = ret_string.split("e")
    # remove leading "+" and strip leading zeros
    sigb = str(int(sigb))
    return siga + ' x 10^{' + sigb + '}'

def peakleft(df,est_peak_left):
    min_difference = abs(est_peak_left - df.energy[0]) 
    for value in df.energy: 
        difference = abs(est_peak_left - value) 
        if difference < min_difference: 
            min_difference = difference 
            peak_left = value 
    return(peak_left)

def peakright(df,est_peak_right):
    min_difference = abs(est_peak_right - df.energy[0]) 
    for value in df.energy: 
        difference = abs(est_peak_right - value) 
        if difference < min_difference: 
            min_difference = difference 
            peak_right = value 
    return(peak_right)

def peakleftwin(df,est_peak_left):
    peak_left = peakleft(df, est_peak_left)
    peak_left_win = peak_left -5
    return(peak_left_win)

def peakrightwin(df,est_peak_right):
    peak_right = peakright(df, est_peak_right)
    peak_right_win = peak_right +5
    return(peak_right_win)

#interested region of a photopeak
def interested_region(df,est_peak_left,est_peak_right,variable):
        peak_left = peakleft(df, est_peak_left)
        peak_right = peakright(df, est_peak_right)
        interested_energy_range = df.loc[(df["energy"] >= peak_left) & (df["energy"] <= peak_right), f"{variable}"]
        return(interested_energy_range)

#def interested_intensity(df,est_peak_left,est_peak_right):
#        peak_left = peakleft(est_peak_left)
#        peak_right = peakright(est_peak_right)
#        interested_intensity_range = df.loc[(df["energy"] >= peak_left) & (df["energy"] <= peak_right), "intensity"]
#        return(interested_intensity_range)

#background fx
def background(df,est_peak_left,est_peak_right):
    peak_left = peakleft(df, est_peak_left)
    peak_right = peakright(df, est_peak_right)
    #intensity_left = int(float(df.loc[(df['energy'] == peak_left), 'intensity'].to_string(index=False)))
    #intensity_right = int(float(df.loc[(df['energy'] == peak_right), 'intensity'].to_string(index=False)))
    #baseline = sum([intensity_left, intensity_right])/len([intensity_left, intensity_right])
    low_mean = df.loc[(df["energy"] >= peak_left - 2) & (df["energy"] <= peak_left), "intensity"].mean()
    high_mean = df.loc[(df["energy"] >= peak_right) & (df["energy"] <= peak_right + 2), "intensity"].mean()
    bg_val = (low_mean + high_mean) / 2 
    return (bg_val)

#peak area fx
def findpeakarea(df,est_peak_left,est_peak_right):
    peak_left = peakleft(df, est_peak_left)
    peak_right = peakright(df, est_peak_right)
    bg_val = background(df,est_peak_left,est_peak_right)
    peak_df = df.loc[(df['energy'] >= peak_left) & (df['energy'] <= peak_right), 'intensity']
    peak_df = peak_df.sub(bg_val)
    peak_sum = peak_df.sum()
    return (peak_sum)

#FWHM fx
def FWHM(df,est_peak_left,est_peak_right):
        energy_interested = interested_region(df, est_peak_left, est_peak_right, "energy")
        intensity_interested = interested_region(df, est_peak_left, est_peak_right, "intensity")
        X = energy_interested.to_numpy()
        Y = intensity_interested.to_numpy()
        bg_val = background(df,est_peak_left,est_peak_right)
        Y_2 = Y - bg_val
        spline = UnivariateSpline(X, Y_2-np.max(Y_2)/2, s=0)
        r1, r2 = spline.roots() # find the roots
        FWHM = r2 - r1
        return(FWHM)


def peakfinder(df,prominence):
    #finding peak and modifying data frame to include peak info
    intensity = df.intensity.to_numpy()
    #print(intensity) 
    peaks = find_peaks(intensity, prominence= prominence)
    #print(peaks)
    peak_index = peaks[0]
    prom = peaks[1]
    #print(prom)
    prominences = prom['prominences']
    #print(prominences)
    my_dict = {
        "peak_index": peak_index,
        "prominences": prominences,
    }
    df_peak = pandas.DataFrame(my_dict)
    df_2 = pandas.merge(df_peak, df, left_on="peak_index", right_index=True, how='outer')

    #this is dataframe that contains only the peak rows
    df_peak = df_2[df_2["prominences"].notnull()]
    #chopping off the first 5 items because they are weird
    peak_energy = df_peak.energy.to_numpy()[5:]

    #finding the list of peak left and peak right
    est_peak_left_arr = []
    est_peak_right_arr = []
    for i in df_peak.energy:
        est_peak_left = int(float((df.loc[(df['energy'] == i), "energy"] - 3).to_string(index=False)))
        est_peak_left_arr.append(est_peak_left)
        est_peak_right = int(float((df.loc[(df['energy'] == i), "energy"] + 3).to_string(index=False)))
        est_peak_right_arr.append(est_peak_right)

    #chopping off the first 5 items because they are weird
    est_peak_left_arr = np.array(est_peak_left_arr)[5:]
    est_peak_right_arr = np.array(est_peak_right_arr)[5:]

    identified_peaks = list(zip(est_peak_left_arr,est_peak_right_arr))
    return(peak_energy,identified_peaks)

def gauss(E, a=1, b=1, c=1):
    sigma = (a + b * (E + c * E**2)**0.5) / (2 * (2 * np.log(2))**0.5)
    return np.random.normal(loc=E, scale=sigma)

def broad_spectrum(pulse_height_values, energy_bin_boarders, sum_intensity, a, b, c):
    energy_bin_centers = (energy_bin_boarders[1:] + energy_bin_boarders[:-1]) / 2
    #energy_bin_centers = energy_bin_boarders
    #a, b, c = -0.643e3, 6.794, 0.258e-3
    number_broadening_samples = sum_intensity

    samples = np.random.choice(energy_bin_centers[0:], size=int(number_broadening_samples), p=pulse_height_values[0:]/np.sum(pulse_height_values[0:]))
    broaded_pulse_height_values = gauss(samples, a, b, c)

    broaded_spectrum, _ = np.histogram(broaded_pulse_height_values, bins=energy_bin_boarders)
    renormalized_broadened_spectrum = broaded_spectrum / np.sum(broaded_spectrum) * np.sum(pulse_height_values[0:])
    
    return(renormalized_broadened_spectrum)

#data_extract functions
def read_file(filepath):
    """ very boring utility function to read a file and create an
        list with each entry a single line from the file
        warning: do not use with very large files (Gb +)
    """
    with open(filepath) as f:
        line_data = f.read().splitlines()
    f.close()

    return line_data

def get_counts(line_data):
    """ extracts the counts from the $ spe file"""
    counts = []
    for i, line in enumerate(line_data):
        if line == "$DATA:":
            startpoint = i + 2

            nchannels = line_data[i+1]
            nchannels = nchannels.split()[-1]

            counts = line_data[startpoint:(startpoint+1+int(nchannels))]

    counts = np.array(counts).astype(int)
    return counts

def get_live_time(line_data):
    """ extracts the live time from the $ spe file"""
    for i, line in enumerate(line_data):
        if line == "$MEAS_TIM:":
            live_time = line_data[i+1]
            live_time = live_time.split()[0]
    return float(live_time)
#package
import openmc
import numpy as np
import pandas
import matplotlib.pyplot as plt

#other files
from functions import peak_area_finder, peak_finder, broaden_spectrum

def foil_peakarea(peakfinder_prominence, energy_bins):
    #extract tallies into pandas df
    with openmc.StatePoint('statepoint.20.h5') as sp:
    #sp = openmc.StatePoint(f"statepoint.20.h5")
        tally = sp.get_tally(name=f"{"pulse-height"}")
    intensity = list(tally.get_values(scores=["pulse-height"]).flatten())
    #chopping first entry because it is ridiculous and last entry because it is excessive
    intensity = intensity[1:]
    energy = energy_bins[1:]
    energy_adjusted = energy[1:]

    df_openmc = pandas.DataFrame({'energy': energy_adjusted, 'intensity': intensity})

    #extract peakfinder results (peak energy, peak energy range)
    peak_energy_arr = peak_finder(df_openmc,peakfinder_prominence)[0]
    identified_peak_range = peak_finder(df_openmc,peakfinder_prominence)[1]

    #array of peak's area
    peak_area_arr = []
    for i,j in identified_peak_range:
        peak_area_arr.append(findpeakarea(df_openmc,i,j))

    df_peak_areas = pandas.DataFrame({'peak_energy': peak_energy_arr, 'peak_areas': peak_area_arr})
    return(df_peak_areas)

def foil_specific_peak_finder(peakfinder_prominence, energy_bins, peak_energy):
    df_peak_area = foil_peakarea(peakfinder_prominence, energy_bins)
    df_peak_area['energy_diff'] = abs(df_peak_area['peak_energy'] - peak_energy) #column of differences between interest peak energy and peak energy
    closest_peak = df_peak_area.loc[df_peak_area['energy_diff'].idxmin()] #find the row which has the smallest differences
    peak_area = closest_peak['peak_areas'] #select the peak area of that energy level
    return(peak_area)

def foil_spectrum(energy_bins):
    with openmc.StatePoint('statepoint.20.h5') as sp:
    #extract tallies into pandas df
        sp = openmc.StatePoint(f"statepoint.20.h5")
    tally = sp.get_tally(name=f"{"pulse-height"}")
    intensity = list(tally.get_values(scores=["pulse-height"]).flatten())
    intensity = intensity[1:]
    #chopping first entry because it is ridiculous and last entry because it is excessive
    energy = energy_bins[1:]
    energy_adjusted = energy[1:]

    df_openmc = pandas.DataFrame({'energy': energy_adjusted, 'intensity': intensity})
    plt.semilogy(df_openmc.energy, df_openmc.intensity)
    plt.xlabel('Energy [keV]')
    plt.ylabel('Intensity')
    plt.show()

def foil_spectrum_processed(energy_bins, peakfinder_prominence, fit_a, fit_b, fit_c):
    #extract tallies into pandas df
    sp = openmc.StatePoint(f"statepoint.20.h5")
    tally = sp.get_tally(name=f"{"pulse-height"}")
    intensity = list(tally.get_values(scores=["pulse-height"]).flatten())
    intensity = intensity[1:]
    #chopping first entry because it is ridiculous and last entry because it is excessive
    energy = energy_bins[1:]
    energy_adjusted = energy[1:]

    df_openmc = pandas.DataFrame({'energy': energy_adjusted, 'intensity': intensity})
    renorm_broadened_intensity = broaden_spectrum(df_openmc.intensity.to_numpy(), energy, sum(df_openmc.intensity), fit_a, fit_b, fit_c) #gaussian broadening
    plt.plot(df_openmc.energy, renorm_broadened_intensity)
    plt.xlabel('Energy [keV]')
    plt.ylabel('Intensity')
    #extract peakfinder results (peak energy, peak energy range)
    peak_energy_arr = peak_finder(df_openmc,peakfinder_prominence)[0]

    for i in peak_energy_arr:
        plt.vlines(x=i, color="red", ls =':', label=f"{i}keV", ymin = 0, ymax=1e15)

    plt.show()
import pandas
import numpy as np
import matplotlib.pyplot as plt

from foil_presimulation import foil_interaction, foil_radionuclide, foil_radionuclide_quantity, foil_gamma
from foil_openmc import foil_simulation
from foil_post_analysis import foil_peakarea, foil_spectrum, foil_spectrum_processed, foil_specific_peak_finder

#main function
def ta181_foil_process(beam_current, cooling_time, counting_time, distance_foil, irradiation_time, proton_energy, thickness_foil, peak_energy):
    cxfile = 'ta181cx.csv'
    foil_isotope = 'ta181'
    #tantulum constasnt
    atomic_mass_foil = 180.94788 #amu
    density_foil = 0.01665 #g/mm^3
    #setup
    energy_bins = np.linspace(0,1000,3000) #tally bins energy
    peakfinder_prominence = 100 #for peakfinder
    fit_a = 1.288356672289505 #these are numbers found for the hpge in eu152_0703 data
    fit_b = 0.04077411501380705
    fit_c = -0.0001950029411875738

    beam_interaction = foil_interaction(atomic_mass_foil, beam_current, density_foil, irradiation_time, thickness_foil) #foil irradiation
    radionuclide_list = foil_radionuclide(cxfile) #name list of foil radionuclides
    df_radionuclide = foil_radionuclide_quantity(proton_energy, cxfile, beam_interaction) #data of foil radionuclides
    df_gamma = foil_gamma(df_radionuclide, radionuclide_list, cooling_time, counting_time, foil_isotope) #creating the gamma energy profile

    foil_simulation(df_gamma, distance_foil, energy_bins) #openmc simulation
    
    peak_area = foil_specific_peak_finder(peakfinder_prominence, energy_bins, peak_energy)
    #foil_spectrum(energy_bins) #display simple spectrum
    #foil_spectrum_processed(energy_bins, peakfinder_prominence, fit_a, fit_b, fit_c) #display processed spectrum (gaussian broadening and peak label)

    return(peak_area)

#test
'''
#variable
beam_current = 1e-6 #A
cooling_time = 100 #s
counting_time = 900 #s
distance_foil = 5 #cm
irradiation_time = 1000 #s
proton_energy = 2e7 #eV
thickness_foil = 0.01 #mm
peak_energy = 433

ta181_foil_process(beam_current, cooling_time, counting_time, distance_foil, irradiation_time, proton_energy, thickness_foil, peak_energy)
'''





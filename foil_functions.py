import pandas
import numpy as np
import matplotlib.pyplot as plt

from foil_presimulation import foil_interaction, foil_radionuclide, foil_radionuclide_quantity, foil_gamma, foil_dose
from foil_openmc import foil_simulation
from foil_post_analysis import foil_peakarea, foil_spectrum, foil_spectrum_processed, foil_specific_peak_finder


def foil_process(parameters, proton_energy, peak_energy):
    #insert tantulum info
    cxfile = 'ta181cx.csv'
    foil_isotope = 'ta181'
    atomic_mass_foil = 180.94788 #amu
    density_foil = 0.01665 #g/mm^3
    thickness_foil = 0.01 #mm
    #setup
    energy_bins = np.linspace(0,1000,3000) #tally bins energy
    peakfinder_prominence = 100 #for peakfinder
    fit_a = 1.288356672289505 #these are numbers found for the hpge in eu152_0703 data
    fit_b = 0.04077411501380705
    fit_c = -0.0001950029411875738

    beam_interaction = foil_interaction(atomic_mass_foil, parameters[0], density_foil, parameters[4], thickness_foil) #foil irradiation
    radionuclide_list = foil_radionuclide(cxfile) #name list of foil radionuclides
    df_radionuclide = foil_radionuclide_quantity(proton_energy, cxfile, beam_interaction) #data of foil radionuclides
    df_gamma = foil_gamma(df_radionuclide, radionuclide_list, parameters[1], parameters[2], foil_isotope) #creating the gamma energy profile
    foil_simulation(df_gamma, parameters[3], energy_bins) #openmc simulation

    peak_area = foil_specific_peak_finder(peakfinder_prominence, energy_bins, peak_energy)
    #foil_spectrum(energy_bins) #display simple spectrum
    #foil_spectrum_processed(energy_bins, peakfinder_prominence, fit_a, fit_b, fit_c) #display processed spectrum (gaussian broadening and peak label)

    return (peak_area)

def foil_dose_process(parameters, proton_energy):
    #insert tantulum info
    cxfile = 'ta181cx.csv'
    foil_isotope = 'ta181'
    atomic_mass_foil = 180.94788 #amu
    density_foil = 0.01665 #g/mm^3
    thickness_foil = 0.01 #mm

    beam_interaction = foil_interaction(atomic_mass_foil, parameters[0], density_foil, parameters[4], thickness_foil) #foil irradiation
    radionuclide_list = foil_radionuclide(cxfile) #name list of foil radionuclides
    df_radionuclide = foil_radionuclide_quantity(proton_energy, cxfile, beam_interaction) #data of foil radionuclides
    dose = foil_dose(df_radionuclide, radionuclide_list, parameters[1], parameters[2], foil_isotope)
    dose_uSv = dose * 1000 * 1.60217e-19 * 1e6 #converting from keV to uJ

    return(dose_uSv)

def foil_process_separate(parameters, proton_energy, peak_energy):
    #insert tantulum info
    cxfile = 'ta181cx.csv'
    foil_isotope = 'ta181'
    atomic_mass_foil = 180.94788 #amu
    density_foil = 0.01665 #g/mm^3
    thickness_foil = 0.01 #mm
    #setup
    energy_bins = np.linspace(0,1000,3000) #tally bins energy
    peakfinder_prominence = 100 #for peakfinder
    fit_a = 1.288356672289505 #these are numbers found for the hpge in eu152_0703 data
    fit_b = 0.04077411501380705
    fit_c = -0.0001950029411875738

    beam_interaction = foil_interaction(atomic_mass_foil, parameters[0], density_foil, parameters[4], thickness_foil) #foil irradiation
    radionuclide_list = foil_radionuclide(cxfile) #name list of foil radionuclides
    df_radionuclide = foil_radionuclide_quantity(proton_energy, cxfile, beam_interaction) #data of foil radionuclides
    df_gamma = foil_gamma(df_radionuclide, radionuclide_list, parameters[1], parameters[2], foil_isotope) #creating the gamma energy profile
    #foil_simulation(df_gamma, parameters[3], energy_bins) #openmc simulation

    peak_area = foil_specific_peak_finder(peakfinder_prominence, energy_bins, peak_energy)
    #foil_spectrum(energy_bins) #display simple spectrum
    foil_spectrum_processed(energy_bins, peak_energy, fit_a, fit_b, fit_c) #display processed spectrum (gaussian broadening and peak label)

    return (peak_area)
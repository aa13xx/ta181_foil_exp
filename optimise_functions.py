import pandas
import numpy as np

from foil_presimulation import foil_interaction, foil_radionuclide, foil_radionuclide_quantity
from main_function import ta181_foil_process

def foil_activity(df_radionuclide, radionuclide_list, cooling_time, counting_time,foil_isotope):
    df_gamma = []
    for i in radionuclide_list:
        df = pandas.read_csv(f"{foil_isotope}_decay/{i}.csv") #will need energy (in kev), intensity (in %), and half_life (in s) **, half_life_m and fraction for daughter
        df['decay'] = 0 #create a new df column that is actual intensity, if
        df['decay'] = df['decay'].astype(float)
        N_1_ti = np.exp((-np.log(2)/df["half_life"]) * cooling_time)
        N_1_tf = np.exp((-np.log(2)/df["half_life"]) * (cooling_time + counting_time))
        N_2_ti = np.exp((-np.log(2)/df["half_life_m"]) * cooling_time) - np.exp((-np.log(2)/df["half_life"]) * cooling_time)
        N_2_tf = np.exp((-np.log(2)/df["half_life_m"]) * (cooling_time + counting_time)) - np.exp((-np.log(2)/df["half_life"]) * (cooling_time + counting_time))
        N_1_ti_ = np.exp((-np.log(2)/df["half_life_m"]) * cooling_time)
        N_1_tf_ = np.exp((-np.log(2)/df["half_life_m"]) * (cooling_time + counting_time))
        lambda_1_2_1 = ((np.log(2)/df["half_life_m"])/((np.log(2)/df["half_life"])-(np.log(2)/df["half_life_m"])))
        df.loc[df['fraction_m'] == 0, 'decay'] = float(df_radionuclide[f"{i}"].iloc[0]) * (N_1_ti - N_1_tf) * 0.01 * df['intensity']
        df.loc[(df['fraction_m'] > 0), 'decay'] = float(df_radionuclide[f"{i}"].iloc[0]) * df['fraction_m'] * (N_1_ti_ - N_1_tf_ - (lambda_1_2_1 * (N_2_tf - N_2_ti))) * 0.01 * df['intensity']
        df_gamma.append(df) 
    df_gamma = pandas.concat(df_gamma,ignore_index=True)
    df_gamma['energy_total'] = df_gamma['decay'] * df_gamma["energy"]
    foil_activity_energy = df_gamma['energy_total'].sum()
    
    return(foil_activity_energy)

def ta181_foil_activity(beam_current, cooling_time, counting_time, distance_foil, irradiation_time, proton_energy, thickness_foil):
    cxfile = 'ta181cx.csv'
    foil_isotope = 'ta181'
    #tantulum constasnt
    atomic_mass_foil = 180.94788 #amu
    density_foil = 0.01665 #g/mm^3

    beam_interaction = foil_interaction(atomic_mass_foil, beam_current, density_foil, irradiation_time, thickness_foil) #foil irradiation
    radionuclide_list = foil_radionuclide(cxfile) #name list of foil radionuclides
    df_radionuclide = foil_radionuclide_quantity(proton_energy, cxfile, beam_interaction) #data of foil radionuclides
    ta181_activity_energy = foil_activity(df_radionuclide, radionuclide_list, cooling_time, counting_time, foil_isotope) 

    return(ta181_activity_energy)


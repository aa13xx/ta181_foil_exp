import pandas
import numpy as np

def check_negative_values(dataframe):
    if (dataframe < 0).any().any():
        raise ValueError("The DataFrame contains negative values.")
    else:
        print("No negative values found in the DataFrame.")

def foil_interaction(atomic_mass_foil, beam_current, density_foil, irradiation_time, thickness_foil):
    avogadros = 6.02214e23
    mm_to_barn = 1e22
    e_charge = 1.60217e-19
    number_density_foil = density_foil * avogadros / atomic_mass_foil 
    beam_interaction = (beam_current / e_charge) * number_density_foil * thickness_foil * irradiation_time / mm_to_barn #b^-1
    return(beam_interaction)

def foil_radionuclide_quantity(proton_energy, cxfile, beam_interaction):
    df_cx = pandas.read_csv(cxfile) #cx data at each proton energy 5-30MeV (columns of reaction product nuclides)
    df_cx = df_cx[(df_cx.energy_p == proton_energy)] #select proton energy
    df_radionuclide = df_cx * beam_interaction 
    return(df_radionuclide)

def foil_radionuclide(cxfile):
    df_cx = pandas.read_csv(cxfile) #cx data at each proton energy 5-30MeV (columns of reaction product nuclides)
    radionuclide_list = df_cx.columns.values[1:]
    #radionuclide_list = radionuclide_list[:2] #this is just for testing if we want to focus on one/limited sets of nuclide gamma
    return(radionuclide_list)

def foil_gamma(df_radionuclide, radionuclide_list, cooling_time, counting_time,foil_isotope):
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
    #print(df_gamma)    
    df_gamma = pandas.concat(df_gamma,ignore_index=True)

    try: #check negative values in df_gamma
        check_negative_values(df_gamma)
    except ValueError as e:
        print(e)

    source_activity = df_gamma['decay'].sum()
    #print(source_activity)
    df_gamma['intensity_real'] =  df_gamma['decay'] / source_activity
    df_gamma = df_gamma.groupby('energy', as_index=False).agg({'decay': 'sum','intensity_real': 'sum','intensity': 'sum'})
    #print(df_gamma)
    return(df_gamma)

def foil_dose(df_radionuclide, radionuclide_list, cooling_time, counting_time,foil_isotope):
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
    foil_dose = df_gamma['energy_total'].sum()
    
    return(foil_dose)

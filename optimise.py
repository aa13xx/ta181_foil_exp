from scipy.optimize import minimize
import pandas
import numpy as np

from optimise_functions import ta181_foil_activity
from main_function import ta181_foil_process

#bounds
beam_current_bounds = (1e-9, 10e-6) #A
cooling_time_bounds = (120, 864000) #s
counting_time_bounds = (0, 86400) #s
distance_foil_bounds = (0.01, 50) #cm
irradiation_time_bounds = (30,3600) #s

boundary = (beam_current_bounds, cooling_time_bounds, counting_time_bounds, distance_foil_bounds, irradiation_time_bounds)

#initial values
beam_current = 10e-7
cooling_time = 10000 #s
counting_time = 1000 #s
distance_foil = 5 #cm
irradiation_time = 1000 #s

initial_guess = (beam_current, cooling_time, counting_time, distance_foil, irradiation_time)

#this is the main objective function which takes in values and spit out the peak area of an interested peak energy
def objective(beam_current, cooling_time, counting_time, distance_foil, irradiation_time):
    peak_energy = 100 #interested peak (in keV)
    #fixed
    proton_energy = 2e7 #eV
    thickness_foil = 0.01 #mm
    peak_area = ta181_foil_process(beam_current, cooling_time, counting_time, distance_foil, irradiation_time, proton_energy, thickness_foil, peak_energy)

    return -peak_area
#constraints
def constraints1(beam_current, cooling_time, counting_time, distance_foil, irradiation_time):
    #fixed
    proton_energy = 2e7 #eV
    thickness_foil = 0.01 #mm
    #constraint values
    foil_activity = 1e12 #in keV
    x = ta181_foil_activity(beam_current, cooling_time, counting_time, distance_foil, irradiation_time, proton_energy, thickness_foil) - foil_activity
     
    return(x) #activity shouldnt exceed the preset foil_activity (radiation amount accumulated during gamma collection in keV)

#constraint_1 = {"type": "ineq", "fun": constraints1}

sol = minimize(objective, initial_guess, method = "SLSQP", bounds = boundary) #, constraints = constraint_1#
print(sol)
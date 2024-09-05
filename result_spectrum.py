from scipy.optimize import minimize
import pandas
import numpy as np

from foil_functions import foil_process_separate
#initial values
beam_current = 1e-5 #A
cooling_time = 120 #s
counting_time = 121.9 #s
distance_foil = 5 #cm
irradiation_time = 124.2 #s

initial_guess = [beam_current, cooling_time, counting_time, distance_foil, irradiation_time]

#Spectrum
proton_energy = 2e7 #eV
peak_energy = 343.4 #keV
print(foil_process_separate(initial_guess, proton_energy, peak_energy))
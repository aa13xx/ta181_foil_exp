from scipy.optimize import minimize
import pandas
import numpy as np

from foil_functions import foil_process, foil_dose_process

#bounds
beam_current_bounds = (1e-9, 10e-6) #A
cooling_time_bounds = (120, 864000) #s
counting_time_bounds = (100, 86400) #s
distance_foil_bounds = (0.1, 20) #cm
irradiation_time_bounds = (30,3600) #s

boundary = [beam_current_bounds, cooling_time_bounds, counting_time_bounds, distance_foil_bounds, irradiation_time_bounds]

#initial values
beam_current = 9e-6 #A
cooling_time = 10000 #s
counting_time = 1000 #s
distance_foil = 5 #cm
irradiation_time = 3600 #s

initial_guess = [beam_current, cooling_time, counting_time, distance_foil, irradiation_time]

#this is the main objective function which aims to minimise the total time
def objective(parameters):
    #parameters 
    beam_current = parameters[0]
    cooling_time = parameters[1]
    counting_time = parameters[2]
    distance_foil = parameters[3]
    irradiation_time = parameters[4]
    
    total_time = cooling_time + counting_time + irradiation_time 
    print(f'parameters = {parameters}')
    print(f'total time = {total_time}')
    return(total_time)

#main constraint = peak area should be = peak_area_goal (10,000) motivated by statistical error
def constraint_main(parameters):
    #parameters 
    beam_current = parameters[0]
    cooling_time = parameters[1]
    counting_time = parameters[2]
    distance_foil = parameters[3]
    irradiation_time = parameters[4]
    #fixed
    proton_energy = 2e7 #eV
    #interested peak
    peak_energy = 343.4 #keV
    peak_area_goal = 10000
    
    x = -peak_area_goal + foil_process(parameters, proton_energy, peak_energy)
    print(f"peak_area_diff = {x}")
    return(x)

#dose constraint
def constraint_dose(parameters):
    #parameters
    beam_current = parameters[0]
    cooling_time = parameters[1]
    counting_time = parameters[2]
    distance_foil = parameters[3]
    irradiation_time = parameters[4]
    #fixed
    proton_energy = 2e7 #eV
    #constraint values
    foil_dose_limit = 50 #mSv/hr
    foil_dose_max = counting_time * foil_dose_limit * 1e-3 / 1.60217e-19 / 3600 / 1000 #converting (mSv/hr) to (keV over counting_time)
    #print(f'foil dose max limit = {foil_dose_max}')
    dose = foil_dose_process(parameters, proton_energy)
    print(f"dose = {dose}")
    x = foil_dose_max - dose
    return(x) 

constraint_dose(initial_guess)

cons_main = {"type": "eq", "fun": constraint_main}
cons_dose = {"type": "ineq", "fun": constraint_dose}
cons = [cons_main, cons_dose]

iteration_results = []

def callback(xk):
    iteration_results.append(xk.copy())

sol = minimize(objective, initial_guess, method = "SLSQP", bounds = boundary, constraints = cons, callback=callback, options={'disp': True})
print(f"solution = {sol}")

df = pandas.DataFrame(iteration_results, columns=['parameter', 'totaltime'])
df.to_csv('optimization_iterations.csv', index=False)

#this is to test the process with the initial values
#print(objective(initial_guess))
#print(constraints1(initial_guess))
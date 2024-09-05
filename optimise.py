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
beam_current = 1e-5 #A
cooling_time = 80 #s
counting_time = 80 #s
distance_foil = 5 #cm
irradiation_time = 80 #s

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
    foil_dose_limit = 50 #uSv/hr
    foil_dose_limit = 50 * 3600 #uSv/s
    dose = foil_dose_process(parameters, proton_energy) / counting_time #uSv/s
    x = foil_dose_limit - dose

    print(f"dose = {dose}")
    return(x) 

constraint_dose(initial_guess)

cons_main = {"type": "eq", "fun": constraint_main}
cons_dose = {"type": "ineq", "fun": constraint_dose}
cons = [cons_main, cons_dose]

#Track iteration Results
df_results = pandas.DataFrame(columns=["Beam_Current", "Cooling_Time", "Counting_Time", "Distance", "Irradiation_Time", "Total_Time"])

def callback(xk):
    fval = objective(xk)  # Calculate the objective function value for the current parameters
    global df_results
    new_row = pandas.DataFrame({"Beam_Current": [xk[0]], "Cooling_Time": [xk[1]], "Counting_Time": [xk[2]], "Distance": [xk[3]], "Irradiation_Time": [xk[4]], 'Total_Time': [fval]})
    df_results = pandas.concat([df_results, new_row], ignore_index=True)

sol = minimize(objective, initial_guess, method = "SLSQP", bounds = boundary, constraints = cons, callback=callback, options={'maxiter': 1000, 'disp': True})
print(f"solution = {sol}")

df_results.to_csv('optimization_results.csv', index=False)


#this is to test the process with the initial values
#print(objective(initial_guess))
#print(constraints1(initial_guess))

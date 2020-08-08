import numpy as np
from scipy.integrate import solve_ivp
from index_populations import *
from ADE_parameters import *
from ADE_functions_parameters import *
from ADE_diff_equation import *


"""
                 Model Dynamics                    

 Function of all differential equations
 (input & output: list of pop order, see above)
"""

def f(t, pop):
    # Initialize
    dpop = [0 for i in np.arange(len(index))]

    lambdat = l(t, pop)
    
    ODEName = {'S^NV': dSNV, 'S^U': dSU, 'S^V': dSV, 'S^NI': dSNI, 'S^PI': dSPI, \
             'S^ADE': dSADE, \
             'E^U': dEU, 'E^V': dEV, 'E^NI': dENI, 'E^PI': dEPI, 'E^ADE': dEADE, \
             'P^U': dPU, 'P^V': dPV, 'P^NI': dPNI, 'P^PI': dPPI, 'P^ADE': dPADE, \
             'I^U+': dIUplus, 'I^U-': dIUminus, 'I^V': dIV, 'I^NI': dINI, 'I^PI': dIPI, 'I^ADE': dIADE, \
             'I^IV': dIIV, 'I^Itilde': dIItilde, 'I^Istar': dIIstar, \
             'L^U+': dLUplus, 'L^U-': dLUminus, 'L^V': dLV, 'L^NI': dLNI, 'L^PI': dLPI, 'L^ADE': dLADE, \
             'L^IV': dLIV, 'L^Itilde': dLItilde, 'L^Istar': dLIstar, 'L^LV': dLLV, \
             'L^Ltilde': dLLtilde, 'L^Istar': dLLstar,\
             'D': D, 'R': R}
   
    for i in compartments:
        notation = i[0]
        if i[1] == 0 :
            dpop[index[notation]] = ODEName[notation](t, pop, lambdat)
        else:
            for j in upperscripts[:i[1]]:
                exception = ((i[0] in ['I','L']) and (j in ['U'])) or (not(i[0] in ['S']) and (j in ['NV']))
                if exception:
                    continue
                notation = i[0] + '^' + j
                if i[2] == 0 :
                    dpop[index[notation]] = ODEName[notation](t, pop, lambdat)
                else:
                    for k in list(range(1,i[2]+1)) :
                        dpop[index[notation + '_' + str(k)]] = ODEName[notation](t, k, pop, lambdat)
    return dpop


"""
                  Solve ODEs                      
"""

pop0 = [0 for i in np.arange(len(index))] # initialisation

# Initial Conditions of the population

pop0[index['S^U']] = 
pop0[index['S^NV']] = 

# Simulation timelapse : the nuumber of days to simulate

days = 

# Solve

sol_pop = solve_ivp(f, [0,days], pop0, method="RK45")


"""
                Save results                     
"""

np.savetxt("ADE_solution_" + name + ".txt", soln.y)



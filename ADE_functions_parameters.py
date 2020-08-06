from ADE_parameters import *
from index_populations import *

# Total number of individuals in the different divisions of population

def popsum(compartment, upperscript, Nerlangs, t, pop):
    x = 0
    for k in range(1, Nerlangs+1):
        x = x + pop[index[compartment + '^' + upperscript + '_' + str(k)]]
    return x
    
# Total number of individuals that can effectively infect susceptible individuals

def popeff(compartment, upperscript, Nerlangs, t, pop):
    return popsum(compartment, upperscript, Nerlangs, t, pop) - popiso(compartment, upperscript, Nerlangs, t, pop) \
        - phome * pophome(compartment, upperscript, Nerlangs, t, pop)
        
# Thus the numbers of symptomatic infections in the fully contagious states and the number of individuals in the late infectious states

def popsick(compartment, upperscript, Nerlangs, t, pop):
    return f_parameter(subscript = 'sick', upperscript = upperscript) * popsum(compartment, upperscript, Nerlangs, t, pop)

# Total number of individuals isolated in general quarantine wards

def Q(t, pop) :
    Q1 = (fUplus_sick * fiso + (1 - fUplus_sick)) * ( popsum('I','U+', NI, t, pop) + popsum('L','U+', NL, t, pop) ) + \
        fiso * (popsick('I', 'V', NI, t, pop) + popsick('L', 'N', NL, t, pop)) + fiso * \ 
        (popsick('I', 'NI', NI, t, pop) + popsick('L', 'NI', NL, t, pop)) 
    Q2 = fiso * (popsick('I', 'PI', NI, t, pop) + popsick('L', 'PI', NL, t, pop)) + \
        fiso * (popsick('I', 'ADE', NI, t, pop) + popsick('L', 'ADE', NL, t, pop)) + \
        fiso * (popsick('I', 'Istar', NI, t, pop) + popsick('L', 'Istar', NL, t, pop)) + \
        fiso * popsick('L', 'Lstar', NL, t, pop)
    return Q1 + Q2

# Total number of individuals isolated in general quarantine wards

def popiso(compartment, upperscript, Nerlangs, t, pop):
    if upperscript == 'U+':
        if (tiso1 <= t <= tiso2) and (Q(t, pop) <= Qmax):
            return (fUplus * fiso + (1 - fUplus_sick))* popsum(compartment, upperscript, Nerlangs, t, pop)
        elif (tiso1 <= t <= tiso2) and (Q(t, pop) > Qmax):
            return (fUplus * fiso + (1 - fUplus_sick))* popsum(compartment, upperscript, Nerlangs, t, pop) * (Qmax / Q(t, pop))
        else 
            return 0
    elif upperscript in ['NI', 'PI', 'ADE', 'Istar']:
        if (tiso1 <= t <= tiso2) and (Q(t, pop) <= Qmax):
            return fiso * popsick(compartment, upperscript, Nerlangs, t, pop)
        elif (tiso1 <= t <= tiso2) and (Q(t, pop) > Qmax):
            return fiso * popsick(compartment, upperscript, Nerlangs, t, pop) * (Qmax / Q(t, pop))
        else 
            return 0
    elif upperscript == 'Lstar' and compartment == 'L':
        if (tiso1 <= t <= tiso2) and (Q(t, pop) <= Qmax):
            return fiso * popsick(compartment, upperscript, Nerlangs, t, pop)
        elif (tiso1 <= t <= tiso2) and (Q(t, pop) > Qmax):
            return fiso * popsick(compartment, upperscript, Nerlangs, t, pop) * (Qmax / Q(t, pop))
        else 
            return 0
    else :
        print("Error: enter a correct values for the upperscript arguments")

# Total number of individuals  isolated at home

def pophome(compartment, upperscript, Nerlangs, t, pop):
    if upperscript == 'U+':
        if (tiso1 <= t <= tiso2) and (Q(t, pop) > Qmax):
            return (fUplus * fiso + (1 - fUplus_sick))* popsum(compartment, upperscript, Nerlangs, t, pop) * (1 - Qmax / Q(t, pop))
        else 
            return 0
    elif upperscript in ['NI', 'PI', 'ADE', 'Istar']:
        if (tiso1 <= t <= tiso2) and (Q(t, pop) > Qmax):
            return popsick(compartment, upperscript, Nerlangs, t, pop) * (1 - fiso * Qmax / Q(t, pop))
        else 
            return 0
    elif upperscript == 'Lstar' and compartment == 'L':
        if (tiso1 <= t <= tiso2) and (Q(t, pop) > Qmax):
            return popsick(compartment, upperscript, Nerlangs, t, pop) * (1 - fiso * Qmax / Q(t, pop))
        else 
            return 0
    else :
        print("Error: enter a correct values for the upperscript arguments")

    
    
    
    
    
    

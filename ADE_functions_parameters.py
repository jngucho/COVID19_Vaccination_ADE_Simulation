from ADE_parameters import *
from index_populations import *

# reproduction number with seasonal fluctuations and effects of general distancing

def R0(t):
    x = R0bar * ( 1 + a * np.cos(2 * np.pi * (t - tROmax) / 365) ) 
    #* (1 - pGen(t))
    return x

# Contact rates in different stages

# denominator

cD = cP * DP + cI * DI + cL * DL

# Contact rate at the prodromal stage P :

def betaP(t):
    return cP * R0(t) / cD

# Contact rate at the prodromal stage I :

def betaI(t):
    return cI * R0(t) / cD

# Contact rate at the prodromal stage L :

def betaL(t):
    return cL * R0(t) / cD

# Contact reduction parameter.

def pCont(t):  
    if tiso1 <= t <= tiso2:
        return pCont
    else:
        return 0

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

# Force of infection is defined by

def l(t, pop):
    lambdaP , lambdaI , lambdaL = 0, 0, 0 
    
    for i in ['U', 'V', 'NI', 'PI', 'ADE']:
        lambdaP = lambdaP + popsum('P', i, NP, t, pop)
    
    for i in ['U-', 'IV', 'Itilde']:
        lambdaI = lambdaI + popsum('I', i, NI, t, pop)
    for i in ['V', 'NI', 'PI', 'ADE', 'Istar']:
        lambdaI = lambdaI + popeff('I', i, NI, t, pop)
        
    for i in ['U-', 'IV', 'Itilde', 'LV', 'Ltilde']:
        lambdaL = lambdaL + popsum('L', i, NL, t, pop)
    for i in ['V', 'NI', 'PI', 'ADE', 'Istar', 'Lstar']:
        lambdaL = lambdaL + popeff('L', i, NL, t, pop)
    
    return (betaP(t) * lambdaP + betaI(t) * lambdaI + betaL(t) * lambdaL) * (1 - pCont(t))
       
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

    
    
    
    
    
    

"""
            Differential equations
                            
Functions with single ODE for certain departments,
"""

from ADE_parameters import *
from index_populations import *
from ADE_functions_parameters import *

# The dynamic of susceptible individuals is described by the following differential equations

def dSNV(t, pop):
    x = - (l(t, pop) + lext) * pop[index['S^NV']] / N
    return x

def dSU(t, pop):
    x = - (l(t, pop) + lext) * pop[index['S^U']] / N - nu * pop[index['S^U']]
    return x
    
def dSV(t, pop):
    x = nu * pop[index['S^U']] - (l(t, pop) + lext) * pop[index['S^V']]/N - alpha \
        *pop[index['S^V']]
    return x

def dSNI(t, pop):
    x = alpha * fNI_S * pop[index['S^V']] - (l(t, pop) + lext) * pop[index['S^NI']] / N 
    return x

def dSPI(t, pop):
    x = alpha * fPI_S * pop[index['S^V']] - (l(t, pop) + lext) * pop[index['S^PI']] / N 
    return x

def dSADE(t, pop):
    x = alpha * fADE_S * pop[index['S^V']] - (l(t, pop) + lext) * pop[index['S^ADE']] / N 
    return x
    

# The dynamic of latent individuals is described by the following differential equations

def dEU(t, k, pop):
    if k == 1:
        x = (l(t, pop) + lext) * pop[index['S^U']] / N - nu * pop[index['E^U_'+str(k)]] \
            - epsilon * pop[index['E^U_'+str(k)]]
        return x
    else:
        x = epsilon * pop[index['E^U_'+str(k-1)]] - epsilon * pop[index['E^U_'+str(k)]] \
            - nu * pop[index['E^U_'+str(k)]]
        return x

def dEV(t, k, pop):
    if k == 1:
        x = (l(t, pop) + lext) * pop[index['S^V']] / N + nu * pop[index['E^U_'+str(k)]] \
            - epsilon * pop[index['E^U_'+str(k)]] - alpha * pop[index['E^V_'+str(k)]]
        return x
    else:
        x = epsilon * pop[index['E^V_'+str(k-1)]] + nu * pop[index['E^U_'+str(k)]] \
            - epsilon * pop[index['E^U_'+str(k)]] - alpha * pop[index['E^V_'+str(k)]]
        return x

def dENI(t, k, pop):
    if k == 1:
        x = (l(t, pop) + lext) * (pop[index['S^NV']] + pop[index['S^NI']]) / N  \
            + alpha * fNI_E * pop[index['E^V_'+str(k)]] - epsilon * pop[index['E^NI_'+str(k)]]
        return x
    else:
        x = epsilon * pop[index['E^NI_'+str(k-1)]] + alpha * fNI_E * \
            pop[index['E^V_'+str(k)]] - epsilon * pop[index['E^NI_'+str(k)]]
        return x

def dEPI(t, k, pop):
    if k == 1:
        x = (l(t, pop) + lext) * pop[index['S^PI']] / N  \
            + alpha * fPI_E * pop[index['E^V_'+str(k)]] - epsilon * pop[index['E^PI_'+str(k)]]
        return x
    else:
        x = epsilon * pop[index['E^PI_'+str(k-1)]] + alpha * fPI_E * \
            pop[index['E^V_'+str(k)]] - epsilon * pop[index['E^PI_'+str(k)]]
        return x

def dEADE(t, k, pop):
    if k == 1:
        x = (l(t, pop) + lext) * pop[index['S^ADE']] / N  \
            + alpha * fADE_E * pop[index['E^V_'+str(k)]] - epsilon * pop[index['E^ADE_'+str(k)]]
        return x
    else:
        x = epsilon * pop[index['E^ADE_'+str(k-1)]] + alpha * fADE_E * \
            pop[index['E^V_'+str(k)]] - epsilon * pop[index['E^ADE_'+str(k)]]
        return x

# The dynamic of prodromal individuals is described by the following differential equations

def dPU(t, k, pop):
    if k == 1:
        x = epsilon * pop[index['E^U_'+str(NE)]]  - varphi * pop[index['P^U_'+str(k)]] \
            - nu * pop[index['P^U_'+str(k)]]
        return x
    else:
        x = varphi * pop[index['P^U_'+str(k-1)]] - varphi * pop[index['P^U_'+str(k)]] \
            - nu * pop[index['P^U_'+str(k)]]
        return x

def dPV(t, k, pop):
    if k == 1:
        x = epsilon * pop[index['E^V_'+str(NE)]]  + nu * varphi * pop[index['P^U_'+str(k)]] \
            - varphi * pop[index['P^V_'+str(k)]] - alpha * pop[index['P^V_'+str(k)]]
        return x
    else:
        x = varphi * pop[index['P^V_'+str(k-1)]] + nu * varphi * pop[index['P^U_'+str(k)]] \
            - varphi * pop[index['P^V_'+str(k)]] - alpha * pop[index['P^V_'+str(k)]]
        return x
        
def dPNI(t, k, pop):
    if k == 1:
        x = epsilon * pop[index['E^NI_'+str(NE)]] + alpha * fNI_P * pop[index['P^V_'+str(k)]] \ 
            - varphi * pop[index['P^NI_'+str(k)]] 
        return x
    else:
        x = varphi * pop[index['P^NI_'+str(k-1)]] + alpha * fNI_P * pop[index['P^V_'+str(k)]] \ 
            - varphi * pop[index['P^NI_'+str(k)]] 
        return x   

def dPPI(t, k, pop):
    if k == 1:
        x = epsilon * pop[index['E^PI_'+str(NE)]] + alpha * fPI_P * pop[index['P^V_'+str(k)]] \ 
            - varphi * pop[index['P^PI_'+str(k)]] 
        return x
    else:
        x = varphi * pop[index['P^PI_'+str(k-1)]] + alpha * fPI_P * pop[index['P^V_'+str(k)]] \ 
            - varphi * pop[index['P^PI_'+str(k)]] 
        return x   
        
def dPADE(t, k, pop):
    if k == 1:
        x = epsilon * pop[index['E^ADE_'+str(NE)]] + alpha * fADE_P * pop[index['P^V_'+str(k)]] \ 
            - varphi * pop[index['P^ADE_'+str(k)]] 
        return x
    else:
        x = varphi * pop[index['P^ADE_'+str(k-1)]] + alpha * fADE_P * pop[index['P^V_'+str(k)]] \ 
            - varphi * pop[index['P^ADE_'+str(k)]] 
        return x        
        
# The dynamic of fully infectious individuals is described by the following differential equations

def dIUplus(t, k, pop):
    if k == 1:
        x = varphi * (fsick + (1 - fsick) * fUplus_I) * pop[index['P^U_'+str(NP)]] \ 
            - gamma * pop[index['I^U+_'+str(k)]] 
        return x
    else:
        x = gamma * pop[index['I^U+_'+str(k-1)]] - gamma * pop[index['I^U+_'+str(k)]]
        return x    

def dIUminus(t, k, pop):
    if k == 1:
        x = varphi * (1- fsick - fUplus_I + fsick * fUplus_I) * pop[index['P^U_'+str(NP)]] \ 
            - gamma * pop[index['I^U-_'+str(k)]] - nu * pop[index['I^U-_'+str(k)]]
        return x
    else:
        x = gamma * pop[index['I^U-_'+str(k-1)]] - gamma * pop[index['I^U-_'+str(k)]] \
            - nu * pop[index['I^U-_'+str(k)]]
        return x 
        
def dIV(t, k, pop):
    if k == 1:
        x = varphi * pop[index['P^V_'+str(NP)]] - gamma * pop[index['I^V_'+str(k)]] \
            - alpha * pop[index['I^V_'+str(k)]]
        return x
    else:
        x = gamma * pop[index['I^V_'+str(k-1)]] - gamma * pop[index['I^V_'+str(k)]] \
            - alpha * pop[index['I^V_'+str(k)]]
        return x        

def dINI(t, k, pop):
    if k == 1:
        x = varphi * pop[index['P^NI_'+str(NP)]] + alpha * fNI_I * pop[index['I^V_'+str(k)]] \
            - gamma * pop[index['I^NI_'+str(k)]]
        return x
    else:
        x = gamma * pop[index['I^NI_'+str(k-1)]] + alpha * fNI_I * pop[index['I^V_'+str(k)]] \
            - gamma * pop[index['I^NI_'+str(k)]]
        return x         

def dIPI(t, k, pop):
    if k == 1:
        x = varphi * pop[index['P^PI_'+str(NP)]] + alpha * fPI_I * pop[index['I^V_'+str(k)]] \
            - gamma * pop[index['I^NI_'+str(k)]]
        return x
    else:
        x = gamma * pop[index['I^PI_'+str(k-1)]] + alpha * fPI_I * pop[index['I^V_'+str(k)]] \
            - gamma * pop[index['I^PI_'+str(k)]]
        return x  
        
def dIADE(t, k, pop):
    if k == 1:
        x = varphi * pop[index['P^ADE_'+str(NP)]] + alpha * fADE_I * pop[index['I^V_'+str(k)]] \
            - gamma * pop[index['I^ADE_'+str(k)]]
        return x
    else:
        x = gamma * pop[index['I^ADE_'+str(k-1)]] + alpha * fADE_I * pop[index['I^V_'+str(k)]] \
            - gamma * pop[index['I^ADE_'+str(k)]]
        return x          

def dIIV(t, k, pop):
    if k == 1:
        x = nu * pop[index['I^U-_'+str(k)]] - gamma * pop[index['I^IV_'+str(NP)]] \
            - alpha * pop[index['I^IV_'+str(k)]] 
        return x
    else:
        x = gamma * pop[index['I^IV_'+str(k-1)]] + nu * pop[index['I^U-_'+str(k)]] \
        - gamma * pop[index['I^IV_'+str(NP)]] - alpha * pop[index['I^IV_'+str(k)]]
        return x 
        
def dIItilde(t, k, pop):
    if k == 1:
        x = alpha * fItilde_I * pop[index['I^IV_'+str(k)]] - gamma * pop[index['I^Itilde_'+str(k)]]
        return x
    else:
        x = gamma * pop[index['I^Itilde_'+str(k-1)]] + alpha * fItilde_I * \
            pop[index['I^IV_'+str(k)]] - gamma * pop[index['I^Itilde_'+str(k)]]
        return x    
        
def dIIstar(t, k, pop):
    if k == 1:
        x = alpha * (1 - fItilde_I) * pop[index['I^IV_'+str(k)]] - gamma * \ 
            pop[index['I^Istar_'+str(k)]]
        return x
    else:
        x = gamma * pop[index['I^Istar_'+str(k-1)]] + alpha * (1 - fItilde_I) * \
            pop[index['I^IV_'+str(k)]] - gamma * pop[index['I^Istar_'+str(k)]]
        return x              

# The dynamic of late infectious individuals is described by the following differential equations
        
def dLUplus(t, k, pop):
    if k == 1:
        x = gamma * pop[index['I^U+_'+str(NI)]] - delta * pop[index['L^U+_'+str(k)]] 
        return x
    else:
        x = delta * pop[index['L^U+_'+str(k-1)]] - delta * pop[index['L^U+_'+str(k)]]
        return x   
        
def dLUminus(t, k, pop):
    if k == 1:
        x = gamma * pop[index['I^U-_'+str(NI)]] - delta * pop[index['L^U-_'+str(k)]] \
            - nu * pop[index['L^U-_'+str(k)]]
        return x
    else:
        x = delta * pop[index['L^U-_'+str(k-1)]] - delta * pop[index['L^U-_'+str(k)]] \
            - nu * pop[index['L^U-_'+str(k)]]
        return x         
        
def dLV(t, k, pop):
    if k == 1:
        x = gamma * pop[index['I^V_'+str(NI)]] - delta * pop[index['L^V_'+str(k)]] \
            - alpha * pop[index['L^V_'+str(k)]]
        return x
    else:
        x = delta * pop[index['L^V_'+str(k-1)]] - delta * pop[index['L^V_'+str(k)]] \
            - alpha * pop[index['L^V_'+str(k)]]
        return x          

def dLNI(t, k, pop):
    if k == 1:
        x = gamma * pop[index['I^NI_'+str(NI)]] + alpha * fNI_L * pop[index['L^V_'+str(k)]] \
            - delta * pop[index['L^NI_'+str(k)]]
        return x
    else:
        x = delta * pop[index['L^NI_'+str(k-1)]] + alpha * fNI_L * pop[index['L^V_'+str(k)]] \
            - delta * pop[index['L^NI_'+str(k)]]
        return x
        
def dLPI(t, k, pop):
    if k == 1:
        x = gamma * pop[index['I^PI_'+str(NI)]] + alpha * fPI_L * pop[index['L^V_'+str(k)]] \
            - delta * pop[index['L^PI_'+str(k)]]
        return x
    else:
        x = delta * pop[index['L^PI_'+str(k-1)]] + alpha * fPI_L * pop[index['L^V_'+str(k)]] \
            - delta * pop[index['L^PI_'+str(k)]]
        return x   
        
def dLADE(t, k, pop):
    if k == 1:
        x = gamma * pop[index['I^ADE_'+str(NI)]] + alpha * fADE_L * pop[index['L^V_'+str(k)]] \
            - delta * pop[index['L^ADE_'+str(k)]]
        return x
    else:
        x = delta * pop[index['L^ADE_'+str(k-1)]] + alpha * fADE_L * pop[index['L^V_'+str(k)]] \
            - delta * pop[index['L^ADE_'+str(k)]]
        return x              
        
def dLIV(t, k, pop):
    if k == 1:
        x = gamma * pop[index['I^IV_'+str(NI)]] - delta * pop[index['L^IV_'+str(k)]] \
            - alpha * pop[index['L^IV_'+str(k)]]
        return x
    else:
        x = delta * pop[index['L^IV_'+str(k-1)]]  - delta * pop[index['L^IV_'+str(k)]] \
            - alpha * pop[index['L^IV_'+str(k)]]
        return x         
                
def dLItilde(t, k, pop):
    if k == 1:
        x = gamma * pop[index['I^Itilde_'+str(NI)]] + alpha * fItilde_L * pop[index['L^IV_'+str(k)]] \
            - delta * pop[index['L^Itilde_'+str(k)]]
        return x
    else:
        x = delta * pop[index['L^Itilde_'+str(k-1)]] + alpha * fItilde_L * pop[index['L^IV_'+str(k)]] \
            - delta * pop[index['L^Itilde_'+str(k)]] 
        return x  
        
def dLIstar(t, k, pop):
    if k == 1:
        x = gamma * pop[index['I^Istar_'+str(NI)]] + alpha * (1 - fItilde_L) * pop[index['L^IV_'+str(k)]] \
            - delta * pop[index['L^Istar_'+str(k)]]
        return x
    else:
        x = delta * pop[index['L^Istar_'+str(k-1)]] + alpha * (1 - fItilde_L) * pop[index['L^IV_'+str(k)]] \
            - delta * pop[index['L^Istar_'+str(k)]]
        return x    
        
def dLLV(t, k, pop):
    if k == 1:
        x = nu * pop[index['L^U-_'+str(k)]] - delta * pop[index['L^LV_'+str(k)]] \
            - alpha * pop[index['L^LV_'+str(k)]]
        return x
    else:
        x = delta * pop[index['L^LV_'+str(k-1)]] + nu * pop[index['L^U-_'+str(k)]] \
            - delta * pop[index['L^LV_'+str(k)]] - alpha * pop[index['L^LV_'+str(k)]] \
        return x         
                
def dLLtilde(t, k, pop):
    if k == 1:
        x = alpha * fLtilde_L * pop[index['L^LV_'+str(k)]] - delta * pop[index['L^Ltilde_'+str(k)]]
        return x
    else:
        x = delta * pop[index['L^Ltilde_'+str(k-1)]] + alpha * fLtilde_L * pop[index['L^LV_'+str(k)]] \
            - delta * pop[index['L^Ltilde_'+str(k)]] 
        return x  
        
def dLIstar(t, k, pop):
    if k == 1:
        x = alpha * (1 - fLtilde_L) * pop[index['L^LV_'+str(k)]] - delta * \
            pop[index['L^Lstar_'+str(k)]]
        return x
    else:
        x = delta * pop[index['L^Lstar_'+str(k-1)]] + alpha * (1 - fLtilde_L) * \
            pop[index['L^LV_'+str(k)]] - delta * pop[index['L^Lstar_'+str(k)]]
        return x                 

# The death toll over time is described by the following equations

def dD(t, pop):
    dD1 = fUplus_sick * fdead * pop[index['L^U+_'+str(NL)]] + fsick * fdead * pop[index['L^V_'+str(NL)]] \
        fsick * fdead * pop[index['L^NI_'+str(NL)]] + fPI_sick * fPI_dead * pop[index['L^PI_'+str(NL)]]
    dD2 = fADE_sick * fADE_dead * pop[index['L^ADE_'+str(NL)]] + fIstar_sick * fIstar_dead * \
        pop[index['L^Istar_'+str(NL)]] + fLstar_sick * fLstar_dead * pop[index['L^Lstar_'+str(NL)]]
    return delta * (dD1 + dD2)

# The dynamic in the recovered population is described by the following equations
def dR(t, pop):
    dR1 = alpha * fR_S * pop[index['S^V']] + alpha * fR_E * popsum('E', 'V', NE, t, pop) + alpha * \
        fR_P * popsum('P', 'V', NP, t, pop) + alpha * fR_I * popsum('I', 'V', NI, t, pop) + alpha * fR_L * popsum('L', 'V', NL, t, pop) \
        + delta * (1 - fUplus_sick * fdead) * pop[index['L^U+_'+str(NL)]]
    dR2 = delta * pop[index['L^U-_'+str(NL)]] + delta * (1 - fsick * fdead) * pop[index['L^V_'+str(NL)]] \ 
        + delta * (1 - fsick * fdead) * pop[index['L^NI_'+str(NL)]] + delta * (1 - fPI_sick * fPI_dead) \
        * pop[index['L^PI_'+str(NL)]] + delta * (1 - fADE_sick * fADE_dead) * pop[index['L^ADE_'+str(NL)]]
    dR3 = delta * pop[index['L^IV_'+str(NL)]] + delta * pop[index['L^Istar_'+str(NL)]] \
        + delta * (1 - fIstar_sick * fIstar_dead) * pop[index['L^Istar_'+str(NL)]] + delta * pop[index['L^LV_'+str(NL)]] \
        + delta * pop[index['L^Ltilde_'+str(NL)]] + delta * (1 - fLstar_sick * fLstar_dead) * pop[index['L^Lstar_'+str(NL)]]
    return dR1 + dR2 + dR3



            

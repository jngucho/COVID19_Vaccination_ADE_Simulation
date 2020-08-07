# Author: Kristina B. Helle
# based on code of the riskgroup model (code extracted from jupyter notebook) with Author = "Nessma A. M. Yousif, Looli Alawam N., Pierre M. Ngougoue N., H. Christian Jr. Tsoungui Obama, Kristan. A. Schneider" Copyright = "Copyright 2020, COVID-19 Risk group model" Credits = "Hochschule Mittweida, University of Applied Sciences" license= "" Version="0.1" Maintainer = "H. Christian Jr. Tsoungui Obama" Email="htsoungl@hs-mitteida.de" Status="Developpement"
# -*- coding: utf-8 -*-


# This is the basic model.
# It imports parameters (l 15) to organise the simulation of different scenarios (each to be saved with own name).
# It saves the result to a textfile (last lines) [and some more things], that can be used for plotting etc.

import numpy as np
from scipy.integrate import solve_ivp
import datetime

# Import parameters
# these are stored in extra files to help keeping track of all scenarios
# they include the parameter "name" to save the result of the scenario
from parameters_new1 import *  # basic scenario (all parameters)
from parameters_new_pdist2 import *  # to be changed for other scenarios (only the parameters that are different)

dummy = -1000000  # technical variable without concrete meaning

# General Functions
# Summarized population: sum of entries of pop from start next length elements
def sumFromLength(pop, start, length):  #
    ind = [pop[i] for i in range(start, start + length)]
    return sum(ind)

# Time tracking (may be omitted)
timestart = datetime.datetime.now()

####################################################
#           Short model description                #
####################################################
# Compartments
#   S susceptible
#   infected:
#   infections have different stages:
#       E early infected
#       P podromal
#       I fully infective
#       L late infections
#   infection can follow 3 parallel versions:
#       _s single infected
#       _l latent multiple infected (in model *)
#       _m multiple infected (in model ~)
#   each stage is subdivided into Erlang stages; same division for all versions
#   D dead
#   R recovered (immune)
# Possible paths:
#   by probability from S to E1, then by rates along the stages
#   from all single infected stages by probability to latent multiple infected, then by rates to multiple infected


#####################################################
#               Model parameters                    #
#####################################################
# Numbers are taken from CovidSim default (covidsim.eu, version 1.1)
# if riskgroups have other values they are indicated as well
# Names are mostly taken from riskgroups; if not it is indicated

### Model settings
# Number of Erlang stages
#! NE = NP = NI = NL = NM = 16  # this is defined in parameters (copied here for explanation - same below).
# number of stages per version
NEPIL = NE + NP + NI + NL

# Population size
#! N = 1000000000

# Duration of simulation
#! days = 365
### Properties of the disease
# Average durations of periods in days
#! DE = 4  # 3.7     # D_E...
#! DP = 1  # 1
#! DI = 5  # 7
#! DL = 5  # 7

## Reproduction
# basic reproduction number
#! R0bar = 3  # ro = 1.2
# amplitude of the seasonal fluctuation
#! a = 0  # 0.5
# day when seasonal fluctuation reaches maximum
#! tROmax = 0  # 300

# reproduction number with seasonal fluctuations and effects of general distancing
def R0(t):
    x = R0bar * (
            1 + a * np.cos(2 * np.pi * (t - tROmax) / 365)
    ) * (1 - pGen(t))
    return x


# Contagiousness at the infectious stages P, I, L
#! cP = 0.5  # shouldn't it be >1? infectivity seems to be highest before symptoms start
#! cI = 1  # 0.98
#! cL = 0.5

# Contact rates in different stages
# denominator
cD = (cP * DP) + DI + (cL * DL)
cD = float(cD)

def betaP(t):
    return cP * R0(t) / cD


def betaI(t):
    return cI * R0(t) / cD


def betaL(t):
    return cL * R0(t) / cD


# Fraction  of  infective  contacts causing super-infection  in  late-infectious state
#! q = 0.5  # [value guessed]

# Fraction of infective contacts with super-infected individuals that cause multiple infections
#! m_l = 0.05  # [value guessed]
#! m_m = 0.1  # [value guessed]

# transition rate from latent super-infected to fully super-infected
alpha = NM / float(DE + DP)

# Probability to get sick
#! fsick_s = 0.82  # 0.1
#! fsick_m = 0.9  # [value guessed]
# Probability of sick person to die
#! fdead_s = 0.003  # 0
#! fdead_m = 0.006  # [value guessed]

### Countermeasures
## Case isolation
# Duration (interval during which measures are in place)
#! tiso1 = 0  # 105  # [value from riskgroups]
#! tiso2 = 365  # 180  # [value from riskgroups]


# Probability that a sick patient is isolated
#! fiso_s = 0.2  # [value from a version of riskgroups]
#! fiso_m = 0.3  # [value guessed]
# Maximum carrying capacity of quarantine ward
#! Qmax = 0
# Relative contact reduction for person in home isolation
#! phome = 0.75  # [value from a version of riskgroups]


## Fractions of respective compartment that is in isolation
# Total number of individuals that are supposed to be isolated
def Q(pop):
    x = fiso_s * fsick_s * (
            + sumFromLength(pop, 1 + NE + NP, NI)  # I_s
            + sumFromLength(pop, 1 + NE + NP + NEPIL, NI)  # I_l
            + sumFromLength(pop, 1 + NE + NP + NI, NL)  # L_s
            + sumFromLength(pop, 1 + NE + NP + NI + NEPIL, NL)  # L_l
    ) + fiso_m * fsick_m * (
                + sumFromLength(pop, 1 + NE + NP + 2 * NEPIL, NI)  # I_m
                + sumFromLength(pop, 1 + NE + NP + NI + 2 * NEPIL, NL)  # L_m
        )
    return x

# Isolation in quarantine ward
# fraction of compartments I, L that are in isolation ward
def factorIso(t, mode, pop):
   Qt = Q(pop)
   y = 0                    # if isolation is not in place, the fraction of isolated is 0
   if tiso1 <= t <= tiso2:  # during this time interval, isolation is in place
        if mode != "m":     # factor of sickness and isolation depends on if single or multiple infection
            y = fsick_s * fiso_s
        else:
            y = fsick_m * fiso_m
        if Qt > Qmax:       # if isolation wards are full, only those who find a place are isolated there
            y = y * Qmax / Qt
   return y


# Home isolation (i.e., their effect)
# fraction of compartments I, L that are in home isolation - actually, their effective number on contacts
def factorHome(t, mode, pop):
    Qt = Q(pop)
    y = 0
    if tiso1 <= t <= tiso2 and Qt > Qmax:
        if mode != "m":
            y = fsick_s * fiso_s
        else:
            y = fsick_m * fiso_m
        y = y * (1 - Qmax / Qt)
    return y


# Not isolated, effectively spreading the disease
# fraction of compartments I, L
def factorEff(t, mode, pop):
    y = 1 - factorIso(t, mode, pop) - phome * factorHome(t, mode, pop)
    return y


## General distancing
# Duration interval
#! tdist1 = 1
#! tdist2 = 165
# Prevented fraction of contacts because of general social-distancing measures
#! pdist = 0.75
# resulting contact reduction

def pGen(t):  # in riskgroups: pcontGe
    if tdist1 <= t <= tdist2:
        return pdist
    else:
        return 0


## Infections from outside per day (lambda_Ext)
#! LExt_s = 1  # LExt
#! LExt_m = 0  # [value guessed]


## Force of infection (lambda)
# [normalized by population size]

# single infection
def lambda_s(t, pop):
    effFactor_s = factorEff(t, "s", pop)
    effFactor_l = factorEff(t, "l", pop)
    effFactor_m = factorEff(t, "m", pop)
    x = LExt_s \
        + betaP(t) * (
                + sumFromLength(pop, 1 + NE, NP)  # P_s sum
                + sumFromLength(pop, 1 + NE + NEPIL, NP) * (1 - m_l)  # P_l sum
                + sumFromLength(pop, 1 + NE + 2 * NEPIL, NP) * (1 - m_m)  # P_m sum
        ) + betaI(t) * (
                + sumFromLength(pop, 1 + NE + NP, NI) * effFactor_s  # I_s^Eff sum
                + sumFromLength(pop, 1 + NE + NP + NEPIL, NI) * effFactor_l * (1 - m_l)  # I_l^Eff sum
                + sumFromLength(pop, 1 + NE + NP + 2 * NEPIL, NI) * effFactor_m * (1 - m_m)  # I_m^Eff sum
        ) + betaL(t) * (
                + sumFromLength(pop, 1 + NE + NP + NI, NL) * effFactor_s  # L_s^Eff sum
                + sumFromLength(pop, 1 + NE + NP + NI + NEPIL, NL) * effFactor_l * (1 - m_l)  # L_l^Eff sum
                + sumFromLength(pop, 1 + NE + NP + NI + 2 * NEPIL, NL) * effFactor_m * (1 - m_m))  # L_m^Eff sum)
    return x / N


# multiple infections
def lambda_m(t, pop):
    effFactor_l = factorEff(t, "l", pop)
    effFactor_m = factorEff(t, "m", pop)
    x = LExt_m \
        + betaP(t) * (
                + sumFromLength(pop, 1 + NE + NEPIL, NP) * m_l  # P_m sum
                + sumFromLength(pop, 1 + NE + 2 * NEPIL, NP) * m_m
        ) \
        + betaI(t) * (
                + sumFromLength(pop, 1 + NE + NP + NEPIL, NI) * effFactor_l * m_l  # I_l^Eff sum
                + sumFromLength(pop, 1 + NE + NP + 2 * NEPIL, NI) * effFactor_m * m_m  # I_m^Eff sum
        ) \
        + betaL(t) * (
                + sumFromLength(pop, 1 + NE + NP + NI + NEPIL, NL) * effFactor_l * m_l  # L_l^Eff sum
                + sumFromLength(pop, 1 + NE + NP + NI + 2 * NEPIL, NL) * effFactor_m * m_m  # L_m^Eff sum
        )
    return x / N


# general force of infection
def l(t, pop):
    x = lambda_s(t, pop) + lambda_m(t, pop)
    return x


################################################
#                 Compartments                 #
################################################
# All compartments are modelled by differential equations.
# They are in one list (pop),
# that is updated in each time step,
# so it always hast as entries the number of individuals per compartment at the current time t

# pop (defined below)
# as the number of Erlang stages is variable, the number of entries and which value is where differs; order:
# S(t),
# E1s(t), ... ENEs(t),
# P1s(t), ... PNPs(t),
# I1s(t), ... INIs(t),
# L1s(t), ... LNLs(t),
# E1l(t), ... ENEl(t),
# P1l(t), ... PNPl(t),
# I1l(t), ... INIl(t),
# L1l(t), ... LNLl(t),
# E1m(t), ... ENEm(t),
# P1m(t), ... PNPm(t),
# I1m(t), ... INIm(t),
# L1m(t), ... LNLm(t),
# R(t), D(t)

## Transitions
# Temporal rates (from one Erlang state and stage to the next) are a list of length 1+ NEPIL and order
# X means that there is an extra entry
# thus pop[i] * rate[i] is the number of individuals moving from compartment i to i * 1
# X, E1s -> E1s, ... ENEs -> P1s, ... PNPs -> I1s, ... INIs -> L1s, ... LNLs -> D+R
# which is identical to the following two
# X, E1l -> E1l, ... ENEl -> P1l, ... PNPl -> I1l, ... INIl -> L1l, ... LNLl -> D+R
# X, E1m -> E1m, ... ENEm -> P1m, ... PNPm -> I1m, ... INIm -> L1m, ... LNLm -> D+R
rate = [0] + [NE / DE] * NE + [NP / DP] * NP + [NI / DI] * NI + [NL / DL] * NL


# Super infections (single infected -> latent multiple infected) are a list of length 1 + NEPIL and order
# X, E1s -> E1l, ..., ENEs -> ENEl, P1s -> P1l, ....................LNLs -> LNLl
def superinf(t, pop):
    x = [l(t, pop)] * (1 + NE + NP) \
        + [l(t, pop) * factorEff(t, "s", pop)] * NI \
        + [l(t, pop) * factorEff(t, "s", pop) * q] * NL
    return x


####################################################
#            Differential equations                #
####################################################
# Functions with single ODE for certain departments, to be combined in pop below

def dS(pop, var):
    x = - var[14] * pop[0]
    return x

# Single infected
def ds1(pop, var, superi):
    x = - rate[1] * pop[1] \
        + var[12] * pop[0] \
        - superi[1] * pop[1]
    return x


def ds(i, pop, superi):
    x = rate[(i - 1)] * pop[i - 1] \
        - rate[i] * pop[i] \
        - superi[i] * pop[i]
    return x


# Latent multiple infected
def dl1(pop, superi):
    i = NEPIL + 1
    x = - rate[i - NEPIL] * pop[i] \
        + superi[i - NEPIL] * pop[i - NEPIL] \
        - alpha * pop[i]
    return x


def dl(i, pop, superi):
    x = rate[(i - 1) - NEPIL] * pop[i - 1] \
        - rate[i - NEPIL] * pop[i] \
        + superi[i - NEPIL] * pop[i - NEPIL] \
        - alpha * pop[i]
    return x


# Multiple infected
def dm1(pop, var):
    i = 2 * NEPIL + 1
    x = - rate[i - 2 * NEPIL] * pop[i] \
        + var[13] * pop[0] \
        + alpha * pop[i - NEPIL]
    return x


def dm(i, pop): #
    x = rate[(i - 1) - 2 * NEPIL] * pop[i - 1] \
        - rate[i - 2 * NEPIL] * pop[i] \
        + alpha * pop[i - NEPIL]
    return x


# Recovered, Dead
def dRecovered(pop): #
    x = rate[NEPIL] * (1 - fsick_s * fdead_s) * (pop[NEPIL] + pop[2 * NEPIL]) \
        + rate[NEPIL] * (1 - fsick_m * fdead_m) * pop[3 * NEPIL]
    return x


def dDead(pop): #
    x = rate[NEPIL] * fsick_s * fdead_s * (pop[NEPIL] + pop[2 * NEPIL]) \
        + rate[NEPIL] * fsick_m * fdead_m * pop[3 * NEPIL]
    return x

####################################################
#                Variables                         #
####################################################
# Variables that are not compartments are kept track of explicitly (column)
# for each day (row) - as solving is not done day-wise, this is just a rough estimate
rec = [[dummy for i in np.arange(15 + 1 + NEPIL)] for j in np.arange(days + 1)]

def variables(t, pop):
    x = [
    R0(t),  #0
    betaP(t),  #1
    betaI(t),  #2
    betaL(t),  #3
    Q(pop),  #4
    factorIso(t, "s", pop),  #5
    factorIso(t, "m", pop),  #6
    factorHome(t, "s", pop),  #7
    factorHome(t, "m", pop),  #8
    factorEff(t, "s", pop),  #9
    factorEff(t, "m", pop),  #10
    pGen(t),  #11
    lambda_s(t, pop),  #12
    lambda_m(t, pop), #13
    l(t, pop)  #14
    ]
    return x


####################################################
#                Model Dynamics                    #
####################################################
# Function of all differential equations
# (input & output: list of pop order, see above)

def f(t, pop):
    # Initialize
    out = [0 for i in np.arange(1 + 3 * NEPIL + 2)]

    ## Variables
    var = variables(t, pop)
    superi = superinf(t, pop)

    ## Compartments
    # Susceptible
    out[0] = dS(pop, var)

    # Single infected
    out[1] = ds1(pop, var, superi)
    for i in range(2, 1 + NEPIL):
        out[i] = ds(i, pop, superi)

    # Latent multiple infected
    out[1 + NEPIL] = dl1(pop, superi)
    for i in range(2 + NEPIL, 1 + 2 * NEPIL):
        out[i] = dl(i, pop, superi)

    # Multiple infected
    out[1 + 2 * NEPIL] = dm1(pop, var)
    for i in range(2 + 2 * NEPIL, 1 + 3 * NEPIL):
        out[i] = dm(i, pop)

    # Recovered
    out[1 + 3 * NEPIL] = dRecovered(pop)
    # Dead
    out[1 + 3 * NEPIL + 1] = dDead(pop)

    # record parameters
    rec[int(t)] = var + superi
    print(t)
    return out


###################################################
#                 Solve ODEs                      #
###################################################
# Initial values
pop0 = [0 for i in np.arange(1 + 3 * NEPIL + 2)]

# All individuals are susceptible
pop0[0] = N - E1_0
pop0[1] = E1_0

# Solve
timestartODE = datetime.datetime.now()
soln = solve_ivp(f, [0,days], pop0, method ="RK45",
                 t_eval = np.arange(0, days),
                 dense_output = True)
timeendODE = datetime.datetime.now()

###################################################
#               Save results                      #
###################################################
np.savetxt("superinfection_solution_" + name + ".txt", soln.y)
# np.savetxt("superinfection_time_" + name + ".txt", soln.t)
np.savetxt("superinfection_parameters_" + name + ".txt", rec)

timeend = datetime.datetime.now()

file1 = open("superinfection_timeStamps_" + name + ".txt", "w")
file1.writelines([str(timestart) + "   \n",
                   str(timestartODE) + "   \n",
                   str(timeendODE) + "   \n",
                   str(timeend) + "   "])
file1.close()

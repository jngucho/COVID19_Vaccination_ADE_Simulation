import numpy as np
import matplotlib.pyplot as plt
from parameters_original import dummy
from parameters_test4 import name  # to be changed for new parameters

plt.ion()
plt.rcParams['figure.figsize'] = 8, 3

# Functions to trim the data before plotting
## sum up different compartments; usually all Erlang stages in one stage of the model
def sumFromLength(pop, start, length):  #
    ind = [pop[i] for i in range(start, start + length)]
    return sum(ind)


## for subplot, find the time of the wave (recovered = susceptible)
def findCross(pop, j, k):
    for i in range(0, pop.shape[1]):
        if pop[j][i] <= pop[k][i]:
            return i


## for variables: keep only true entries, delete the ones that have not been changed after initialization
def selectTrueVal(x):
    y = [0 for i in np.arange(len(x))]
    k = 0
    for i in range(0, len(x)):
        if x[i] != dummy:
            y[k] = x[i]
            k += 1
    z = y[0:k]
    return z


# load data
# solution of ODE, model dynamics
pop = np.loadtxt('superinfection_solution_' + name + ".txt")
# variables
var = np.loadtxt('superinfection_parameters_' + name + ".txt")

# determine numbers relevant to select data for plotting
NE = NP = NI = NL = NM = round((len(pop) - 3) / 12)
NEPIL = NE + NP + NI + NL
j = 0
k = 1 + 3 * NEPIL

# compartments; sum up all Erlang stages per stage
popSum = [
    sumFromLength(pop, 1 , NE),  # E_s
    sumFromLength(pop, 1 + NEPIL, NE),  # E_l
    sumFromLength(pop, 1 + 2 * NEPIL, NE),  # E_m
    sumFromLength(pop, 1 + NE, NP),  # P_s
    sumFromLength(pop, 1 + NE + NEPIL, NP),  # P_l
    sumFromLength(pop, 1 + NE + 2 * NEPIL, NP),  # P_m
    sumFromLength(pop, 1 + NE + NP, NI),  # I_s
    sumFromLength(pop, 1 + NE + NP + NEPIL, NI),  # I_l
    sumFromLength(pop, 1 + NE + NP + 2 * NEPIL, NI),  # I_m
    sumFromLength(pop, 1 + NE + NP + NI, NL),  # L_s
    sumFromLength(pop, 1 + NE + NP + NI + NEPIL, NL),  # L_l
    sumFromLength(pop, 1 + NE + NP + NI + 2 * NEPIL, NL)  # L_m
]

# plot compartments
plt.figure()
plt.plot(pop[0], label = "S", color = 'm')
plt.plot(popSum[0], label="Es", color = "y")
plt.plot(popSum[1], label="El", color = "y", linestyle = '--')
plt.plot(popSum[2], label="Em", color = "y", linestyle = ':')
plt.plot(popSum[3], label="Ps", color="g")
plt.plot(popSum[4], label="Pl", color="g", linestyle = '--')
plt.plot(popSum[5], label="Pm", color="g", linestyle = ':')
plt.plot(popSum[6], label="Is", color="c")
plt.plot(popSum[7], label="Il", color="c", linestyle = '--')
plt.plot(popSum[8], label="Im", color="c", linestyle = ':')
plt.plot(popSum[9], label="Ls", color="b")
plt.plot(popSum[10], label="Ll", color="b", linestyle = '--')
plt.plot(popSum[11], label="Lm", color="b", linestyle = ':')
plt.plot(pop[1 + 3 * NEPIL], label="R", color="grey")
plt.plot(pop[1 + 3 * NEPIL + 1], label="D", color="k")
plt.legend(loc= "right")
plt.xlabel ("Days")
plt.ylabel("Individuals")
plt.savefig("Superinfection_" + name + ".png")

# plot subplot of the wave
mid = findCross(pop, j, k)
# mid = 60
plt.ylim([0, pop[0, 0]/5])
plt.xlim([mid - 30, mid + 30])
plt.savefig("Superinfection_subplot_" + name + ".png")

# plot variables
plt.subplot(3, 2, 1)
plt.plot(selectTrueVal(var[:, 0]), label="R0", color="y")
plt.legend(loc= "right")

plt.subplot(3, 2, 2)
plt.plot(selectTrueVal(var[:, 1]), label="betaP", color="k")
plt.plot(selectTrueVal(var[:, 2]), label="betaI", color="k", linestyle = '--')
plt.plot(selectTrueVal(var[:, 3]), label="betaL", color="k", linestyle = ':')
plt.legend(loc= "right")

plt.subplot(3, 2, 3)
plt.plot(selectTrueVal(var[:, 4]), label="Q", color="y")
plt.legend(loc= "right")

plt.subplot(3, 2, 4)
plt.plot(selectTrueVal(var[:, 5]), label="f_iso_s", color="g")
plt.plot(selectTrueVal(var[:, 6]), label="f_iso_m", color="g", linestyle = '--')
plt.plot(selectTrueVal(var[:, 7]), label="f_home_s", color="c")
plt.plot(selectTrueVal(var[:, 8]), label="f_home_m", color="c", linestyle = '--')
plt.plot(selectTrueVal(var[:, 9]), label="f_eff_s", color="b")
plt.plot(selectTrueVal(var[:, 10]), label="f_eff_m", color="b", linestyle = '--')
plt.legend(loc= "right")

plt.subplot(3, 2, 5)
plt.plot(selectTrueVal(var[:, 11]), label="p_ges", color="m", linestyle = ':')
plt.legend(loc= "right")

plt.subplot(3, 2, 6)
plt.plot(selectTrueVal(var[:, 12]), label="lambda_s", color="r")
plt.plot(selectTrueVal(var[:, 13]), label="lambda_m", color="g", linestyle = '--')
plt.plot(selectTrueVal(var[:, 14]), label="lambda", color="b", linestyle = ':')
plt.legend(loc= "right")

plt.savefig("Superinfection_variables_" + name + ".png")






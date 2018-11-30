# This program calculates the total energy and magnetization
# for a 1D Ising model with N dipoles
# Author: Nico Grisouard, University of Toronto
# Date: 20 November 2018

# import modules
import numpy as np
import matplotlib.pyplot as plt
from random import random, randrange


def energyfunction(dipoles):
    """ Function to calculate energy
    IN: dipoles, (4, )-shaped numpy array of ints
    OUT: float of energy"""
    energy = -np.sum(dipoles[0:-2]*dipoles[1:-1])
    return energy


#
def acceptance(En, Eo, kT):
    """ Function for acceptance probability
    IN:
    En [float] the new energy
    Eo [float] the old energy
    kT [float] kB*T
    OUT: accepted [bool]
    """
    p = np.exp(-(En - Eo)/kT)  # Boltzmann factor
    if En-Eo <= 0 or random() < p:
        accepted = True
    else:  # rejected
        accepted = False
    return accepted


def M(dipoles): return np.sum(dipoles)  # total magnetization


# define constants
kB = 1.0
T = 1.0
J = 1.0
num_dipoles = 100
N = 10000  # number of flips

# generate array of dipoles
dipoles = np.ones(num_dipoles, int)
energy = []
magnet = []

E = J*energyfunction(dipoles)
energy.append(E)
magnet.append(M(dipoles))

for i in range(N):
    picked = randrange(num_dipoles)
    dipoles[picked] *= -1  # We flip
    Enew = J*energyfunction(dipoles)

    # calculate new energy depending on probability
    flipd = acceptance(Enew, E, kB*T)  # this is the next old value
    if flipd:
        E = Enew
    else:
        dipoles[picked] *= -1  # we de-flip

    # store energy and magnetization
    energy.append(E)
    magnet.append(M(dipoles))

# plot energy, magnetization
plt.figure()

plt.subplot(211)
plt.plot(energy)
plt.grid()
plt.xlabel("Number of flips")
plt.ylabel("Total energy")

plt.subplot(212)
plt.plot(magnet)
plt.grid()
plt.xlabel("Number of flips")
plt.ylabel("Total magnetization")

plt.tight_layout()
plt.savefig('../images/1D_ising.png')
plt.show()

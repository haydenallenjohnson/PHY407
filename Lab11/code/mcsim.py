# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from random import random,randrange
from math import exp,pi
import numpy as np
import matplotlib.pyplot as plt

T = 10
N = 1000
steps = 250000

# Create a 2D array to store the quantum numbers
n = np.ones([N,3],int)

# Main loop
eplot = []
E = 3*N*pi*pi/2
for k in range(steps):

    # Choose the particle and the move
    i = randrange(N)
    j = randrange(3)
    if random()<0.5:
        dn = 1
        dE = (2*n[i,j]+1)*pi*pi/2
    else:
        dn = -1
        dE = (-2*n[i,j]+1)*pi*pi/2

    # Decide whether to accept the move
    if n[i,j]>1 or dn==1:
        if random()<exp(-dE/T):
            n[i,j] += dn
            E += dE

    eplot.append(E)

# Make the graph
plt.figure(0)
plt.plot(eplot)
plt.title('Energy of the system over time')
plt.ylabel("Energy")
plt.xlabel("Number of steps")
plt.savefig('../images/q1a_energy.png')

#This calculates the energy of each particle, neglecting constant factors
energy_n = n[:,0]**2 + n[:,1]**2 + n[:,2]**2

# This calculates the frequency distribution and creates a plot 
plt.figure(1)
plt.clf()
plt.title('Histogram of normalized energy')
plt.ylabel('Frequency')
plt.xlabel('Normalized energy $e_n = n^2$')
hist_output = plt.hist(energy_n, 50)
plt.savefig('../images/q1a_energy_hist.png')

# This is the frequency distribution
energy_frequency = hist_output[0]
# This is what the x-axis of the plot should look like
# if we plot the energy distribution as a function of n
# the 2nd axis of hist_output contains the boundaries of the array. # Instead, we want their central value, below.
energy_vals = 0.5*(hist_output[1][:-1] + hist_output[1][1:])
n_vals = energy_vals**0.5

# Create the desired plot
plt.figure(2)
plt.clf()
plt.bar(n_vals, energy_frequency, width=0.1)
plt.title('Histogram of values of $n$')
plt.ylabel('Frequency')
plt.xlabel('Value of $n=\sqrt{e_n}$')
plt.savefig('../images/q1a_n_hist.png')

#compute average value of n
n_av = np.sum(energy_frequency*n_vals)/np.sum(energy_frequency)
print('Average value of n:',n_av)

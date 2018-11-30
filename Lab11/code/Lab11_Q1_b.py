#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 13:19:17 2018

@author: Hayden
"""

#import modules
from random import random,randrange
from math import exp,pi
import numpy as np
import matplotlib.pyplot as plt

#define values
T = [10,40,100,400,1200,1600]
N = 1000
equil_cond = 200000

#initialize lists
#eplots = []
n_averages = np.empty(len(T))
E_averages = np.empty(len(T))

#initialize figure
plt.figure()

#iterate over temperature valeus
for l in range(len(T)):
    
    # Create a 2D array to store the quantum numbers
    n = np.ones([N,3],int)

    #initialize main loop
    run = True
    steps = 0
    E_av_old = 0
    
    eplot = []
    E = 3*N*pi*pi/2
    
    #while run condition is true
    while run == True:
        steps += 1
        
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
            if random()<exp(-dE/T[l]):
                n[i,j] += dn
                E += dE

        eplot.append(E)

        #check if equilibrated
        if steps%equil_cond == 0:
            E_av_new = np.mean(eplot[steps-equil_cond:steps])
            
            #if equilibrium reached, store average and break loop
            if E_av_new <= E_av_old:
                run = False
                E_averages[l] = (E_av_new+E_av_old)/2.0
                
            #swap old and new average values
            E_av_old,E_av_new = E_av_new,E_av_old
            
    energy_n = n[:,0]**2 + n[:,1]**2 + n[:,2]**2
    hist_output = plt.hist(energy_n, 50)

    energy_frequency = hist_output[0]
    energy_vals = 0.5*(hist_output[1][:-1] + hist_output[1][1:])
    n_vals = energy_vals**0.5

    n_av = np.sum(energy_frequency*n_vals)/np.sum(energy_frequency)

    n_averages[l] = n_av
    
    plt.plot(eplot, label='$k_BT$='+str(T[l]))
    
plt.grid()
plt.title('Energy as a function of time for system at various temperatures')
plt.xlabel('Number of steps')
plt.ylabel('Energy')
plt.legend(loc='lower right')
plt.xticks([0,0.5e6,1e6,1.5e6,2e6])
plt.savefig('../images/q1b_energies.png')

plt.figure()
plt.plot(T, n_averages)
plt.grid()
plt.title('Average value of $n$ as a function of temperature')
plt.xlabel('$k_BT$')
plt.ylabel('Average value of $n$')
plt.savefig('../images/q1b_av_n.png')

plt.figure()
plt.plot(T,E_averages)
plt.grid()
plt.title('Average value of equilibrium energy as a function of temperature')
plt.xlabel('$k_BT$')
plt.ylabel('Equilibrium energy')
plt.savefig('../images/q1b_av_energy.png')
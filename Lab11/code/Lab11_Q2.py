#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 23:55:11 2018

@author: Hayden
"""

#import modules
import numpy as np
import matplotlib.pyplot as plt

#define function to compute energy of the system
def energy(dipoles):
    vert = np.sum(dipoles[:-1,:]*dipoles[1:,:])
    horiz = np.sum(dipoles[:,:-1]*dipoles[:,1:])
    return -J*(vert + horiz)

# define constants
k_B = 1.0
T = 1.0
J = 1.0
grid_size = 20

#specify number of flibs
N = 200000

#initialize arrays to store magnetization
magnetization = np.empty(N)

#initialize array of dipoles
dipoles = np.random.choice([-1,1], size=(grid_size, grid_size))
E_old = energy(dipoles)

#initialze loop
for k in range(N):
    
    #randomly pick the position of a dipole to flip
    position = np.random.randint(grid_size, size=2)
    i = position[0]
    j = position[1]
    
    #flip the dipole
    dipoles[i,j] *= -1
    
    #calculate the new energy
    E_new = energy(dipoles)
    
    #decide whether or not to accept the move
    prob = np.exp(-(E_new - E_old)/(k_B*T))
    
    if E_new <= E_old or np.random.random() < prob:
        #if move is accepted, update energy value
        E_old,E_new = E_new,E_old
        
    else:
        #if move is rejected, put flipped dipole back how it was
        dipoles[i,j] *= -1
    
    #compute the total magnetization and update the array
    magnetization[k] = np.sum(dipoles)
    
#plot magnetization over time
plt.figure()
plt.plot(magnetization)
plt.grid()
plt.xlabel('Number of steps')
plt.ylabel('Total magnetization')
plt.title('Magnetization of the system over time')
plt.savefig('../images/q2_magnetization.png')
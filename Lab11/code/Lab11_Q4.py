# -*- coding: utf-8 -*-
"""
Uses a simulated annealing method to implement a dimer covering in a 50x50 grid
space and plots the final grid covering.

Created on Wed Nov 28 12:00:00 2018
@author: Pierino Zindel
"""
#import libraries
import numpy as np
import pylab as plt
from random import random, randint, seed

#dimer object
class dimer:
    def __init__(self):    
        """initialize a dimer object"""
        #generate random 1st point
        x1, y1 = randint(1,L), randint(1,L)
        #generate valid random adjacent point
        x2, y2 = self.adj_point(x1,y1)
        while x2 < 1 or x2 > L or y2 < 1 or y2 > L:
            x2, y2 = self.adj_point(x1,y1)
        #assign the points
        self.pt1 = [x1, y1]
        self.pt2 = [x2, y2]
    
    def adj_point(self, x, y):
        """generate valid random adjacent point"""
        #set to the original point
        x2 = x
        y2 = y
        #offset 2nd point accordingly
        n = randint(1,4)
        if n == 1:
            x2 += 1
        elif n == 2:
            x2 -= 1
        elif n == 3:
            y2 += 1
        elif n == 4:
            y2 -= 1
        return x2, y2
        
    def xcoord(self):
        """return the x coordinates of the points in the dimer"""
        return [self.pt1[0], self.pt2[0]]
    
    def ycoord(self):
        """return the y coordinates of the points in the dimer"""
        return [self.pt1[1], self.pt2[1]]


#parameters of the problem
L = 50
steps = 100
Tmax = 10
Tmin = 1e-3
tau = 1e4
Energies = []
times = []

#containder dictionary for all dimers
dimers = {}

#seed for fixed solution
seed(1)

#main loop
t = 0
T = Tmax
while T > Tmin:    
    #compute the 'temperature'
    t += 1
    T = Tmax*np.exp(-t/tau)

    #variables for the current dimer    
    new_dimer = dimer()
    occupied = False
    occupied_by_one = False
    occupation_key = 0
    
    #check how the space is occupied
    for key in dimers.keys():
        fixed = dimers[key]
        #check if any of the spaces are occupied
        if new_dimer.pt1 == fixed.pt1 or new_dimer.pt1 == fixed.pt2 or \
            new_dimer.pt2 == fixed.pt1 or new_dimer.pt2 == fixed.pt2:
                occupied = True
        #check if both spaces occupied by one dimer
        if (new_dimer.pt1 == fixed.pt1 and new_dimer.pt2 == fixed.pt2) or \
            (new_dimer.pt1 == fixed.pt2 and new_dimer.pt2 == fixed.pt1):
                occupied_by_one = True
                occupation_key = key
    
    #assign the grid space accordingly
    if not occupied:
        dimers[t] = new_dimer
    elif occupied_by_one:
        if random() < np.exp(-1/T):
            del dimers[key]
        
    #compute the engery of the system
    E = -len(dimers)
    Energies.append(E)
    times.append(t)

#plot the filled out dimer covering
plt.figure()
plt.title("Dimer covering, Tau=10,000")
plt.ylabel("y position")
plt.xlabel("x position")
plt.xlim(0,L+1)
plt.ylim(0,L+1)
for key in dimers.keys():
    plt.plot(dimers[key].xcoord(), dimers[key].ycoord(), 'o-', markersize=3, linewidth=2, c='k')
plt.grid()
plt.savefig("../images/q4_dimers_t=1e4.png", dpi=500)

#plot the energy over time function
plt.figure()
plt.title("Energy function of dimer covering, Tau=10,000")
plt.ylabel("Energy")
plt.xlabel("Time")
lbl = "$E_{final}=$" + str(Energies[-1])
plt.plot(times, Energies, label=lbl)
plt.grid()
plt.legend()
plt.savefig("../images/q4_energy_t=1e4.png", dpi=500)

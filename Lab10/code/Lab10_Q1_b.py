# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 10:10:20 2018

@author: Pierino
"""
#import libraries
import numpy as np
import pylab as plt
from time import time
from random import randint, seed

#parameters of our problem
L = 101 #grid length
N = 5000 #number of particles


#function that randomly moves the particle and returns the coordinates
def move(i,j):
    #seed and generate a random number
    seed(time())
    n = randint(1,4)
    
    #move the position using the random nunber
    if n == 1:
        i += 1
    elif n == 2:
        i -= 1
    elif n == 3:
        j += 1
    elif n == 4:
        j -= 1
    else:
        raise Exception("invalid move for particle; random value outside of \
                        range")
        
    return i,j


#function to check the boundary and returns turn if the particle becomes stuck
def stuck_check(i,j,particles):
    #check for the boundary
    if i  == 0 or j == 0 or i == (L-1) or j == (L-1):
        return True
    
    #check for surrounding particles
    for particle in particles:
        if i == (particle[0]+1) and j == particle[1]:
            return True
        elif i == (particle[0]-1) and j == particle[1]:
            return True
        elif i == particle[0] and j == (particle[1]+1):
            return True
        elif i == particle[0] and j == (particle[1]-1):
            return True
    
    return False


#plotting function
def plot_traj(x,y):
    plt.figure()
    plt.plot(x, y, '.')
    plt.title("2D Brownian Motion Simulation")
    plt.ylabel("y position")
    plt.xlabel("x position")
    plt.grid()
    plt.show()
    
#container for positions of all stuck particles
particles = []
#container for the trajectories of all particles
trajectories = []

for k in range(N):
    #set starting position
    i = L//2
    j = L//2
    
    #containers for the particles trajectory
    x = [i]
    y = [j]
   
    stuck = False
    #particle_num = len(particles)+1
    
    while not stuck:
        #move the position
        i,j = move(i,j)
        #store the coordinates
        x.append(i)
        y.append(j)
        #check if the particle is done moving
        stuck = stuck_check(i,j,np.array(particles))
        
        #plot the steps
        #plot_traj(x,y)
        
    #store the final coordinates of the particle
    particles.append([i,j])
    #store the full trajectory of the particle
    x = np.array(x)
    y = np.array(y)
    trajectories.append([x,y])
    
    temp = np.array(particles)
    
    plot_traj(temp[:,0],temp[:,1])
    

#plot trajectory of the particle

#plt.savefig("../images/Brownian_moition.png", dpi=500)

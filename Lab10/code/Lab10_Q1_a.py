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
L = 101
N = 5000

#function that randomly moves the particle and returns the coordinates
def move(i,j):
    #seed and generate a random number
    seed(time())
    n = randint(1,4)
    
    #move the position using the random nunber
    if n == 1 and (i+1) < L:
        i += 1
    elif n == 2 and (i-1) >= 0:
        i -= 1
    elif n == 3 and (j+1) < L:
        j += 1
    elif n == 4 and (j-1) >= 0:
        j -= 1
    else:
        i,j = move(i,j) #if updated position is outside of bounds retry move
        
    return i,j

#set starting position
i = L//2
j = L//2

#containers for the particles trajectory
x = [i]
y = [j]

for k in range(N):
    i,j = move(i,j)
    x.append(i)
    y.append(j)
    
#cast containers to arrays
x = np.array(x)
y = np.array(y)

#plot trajectory of the particle
plt.figure()
plt.plot(x, y, label="trajectory")
plt.title("2D Brownian Motion Simulation")
plt.ylabel("y position")
plt.xlabel("x position")
plt.grid()
plt.legend()
#plt.savefig("../images/Brownian_moition.png", dpi=500)

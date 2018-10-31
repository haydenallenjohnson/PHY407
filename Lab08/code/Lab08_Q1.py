#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 31 10:09:32 2018

@author: Hayden
"""

import numpy as np
import pylab as pl

#specify grid size and resolution
height = 8
width = 20
dx = 0.1

#specify target accuracy
target = 1e-6

#create array to hold values
temp = np.zeros((int(height/dx)+1,int(width/dx)+1),float)

#specify boundary conditions of array
temp[-1,:51] = np.linspace(0,5,51)
temp[50:,50] = np.linspace(7,5,31)
temp[50,50:151] = 7
temp[50:,150] = np.linspace(7,5,31)
temp[-1,150:] = np.linspace(5,0,51)
temp[:,-1] = np.linspace(10,0,81)
temp[0,:] = 10
temp[:,0] - np.linspace(10,0,81)

'''
# Main loop
delta = 1.0
while delta>target:

    # Calculate new values of the potential
    for i in range(M+1):
        for j in range(M+1):
            if i==0 or i==M or j==0 or j==M:
                phiprime[i,j] = phi[i,j]
            elif j==M//2 and i>=M//4:
                phiprime[i,j] = phi[i,j]
            else:
                phiprime[i,j] = (phi[i+1,j] + phi[i-1,j] \
                                 + phi[i,j+1] + phi[i,j-1])/4

    # Calculate maximum difference from old values
    delta = max(abs(phi-phiprime))

    # Swap the two arrays around
    phi,phiprime = phiprime,phi

# Make a plot
imshow(phi)
gray()
show()

print(phi[M//4,M//4])
'''
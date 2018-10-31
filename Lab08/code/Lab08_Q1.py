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
w = 0.9

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
temp[:,0] = np.linspace(10,0,81)

#
counter = 0
while counter < 100:
    counter += 1
    # Calculate new values of the potential
    for i in range(temp.shape[0]):
        for j in range(temp.shape[1]):
            if i==0 or i==temp.shape[0]-1 or j==0 or j==temp.shape[1]-1:
                continue
            elif 50<=j<=150 and i>=50:
                continue
            else:
                temp[i,j] = ((1+w)/4)*(temp[i+1,j] + temp[i-1,j] + temp[i,j+1] + temp[i,j-1]) - w*temp[i,j]

    # Calculate maximum difference from old values
#    delta = max(abs(phi-phiprime))
    pl.clf()
    pl.imshow(temp)
    pl.gray()
    pl.pause(0.01)
    pl.show()

# Make a plot



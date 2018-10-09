#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 20 13:53:17 2018

@author: Hayden
"""

#import modules
import numpy as np
import matplotlib.pyplot as plt
from Lab02_Q2_functions import I

#create x and y coordinates of grid (in um)
x = np.arange(-1, 1, 0.01)
y = np.arange(-1, 1, 0.01)

#create empty 2D array to store r values for points on grid
r = np.zeros((len(y),len(x)))
for i in range(len(x)):
    for j in range(len(y)):
        r[i,j] = np.sqrt(x[i]*x[i] + y[j]*y[j])

#compute the intensity for each point on the grid
intensity = I(r, 0.5)

#plot intensity
plt.figure(0)
plt.pcolormesh(x, y, intensity, vmax = 0.01)
plt.suptitle('Plot of light intesity for a diffraction pattern with $\lambda = 500$nm')
plt.xlabel('x ($\mu$m)')
plt.ylabel('y ($\mu$m)')
plt.xlim([-1,1])
plt.ylim([-1,1])
plt.gca().set_aspect('equal', adjustable='box')
plt.savefig('diffraction_intensity.png')

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 12 10:49:58 2018

@author: Hayden
"""
'''
Code for part d) - numerical integration of equations of motion for
mercury under the influence of gravity as described by general 
relativity
'''

#import modules
import numpy as np
import matplotlib.pyplot as plt

#define timestep and interval length
dt = 0.0001
N = 10001

#create time array (for plotting later)
t = np.linspace(0.0, 1.0, N)

#define constants (in units of AU, solar masses, and Earth years)
G = 39.5
M_s = 1
alpha = 0.01

#define empty arrays
x = np.empty(N)
y = np.empty(N)
v_x = np.empty(N)
v_y = np.empty(N)

#set initial positions and velocities
x[0] = 0.47
y[0] = 0.0
v_x[0] = 0.0
v_y[0] = 8.17

#initiate for loop
for i in range(N-1):
    
    #calculate r^2 and r^3 once first to save time
    r_sq = x[i]*x[i] + y[i]*y[i]
    r_cubed = r_sq**(1.5)
    
    #compute GR correction factor
    corr = (1+(alpha/r_sq))
    
    #update velocities
    v_x[i+1] = v_x[i] - corr*G*M_s*x[i]*dt/r_cubed
    v_y[i+1] = v_y[i] - corr*G*M_s*y[i]*dt/r_cubed

    #update positions
    x[i+1] = x[i] + v_x[i+1]*dt
    y[i+1] = y[i] + v_y[i+1]*dt

#calculate angular momenta
L = x*v_y - y*v_x

#plot positions
plt.figure(0)
plt.plot(x,y)
plt.grid()
plt.xlim(-0.5, 0.5)
plt.ylim(-0.5, 0.5)
plt.gca().set_aspect('equal', adjustable='box')
plt.title('Relativistic orbit of Mercury (exaggerated)')
plt.xlabel('x position (AU)')
plt.ylabel('y position (AU)')
plt.savefig('gr_orbit.png')

#plot x velocity vs. time
plt.figure(1)
plt.plot(t, v_x)
plt.grid()
plt.title('x velocity vs. time for relativistic orbit')
plt.ylabel('x velocity (AU/yr)')
plt.xlabel('time (yr)')
plt.savefig('gr_xvel.png')

#plot y velocity vs. time
plt.figure(2)
plt.plot(t, v_y)
plt.grid()
plt.title('y velocity vs. time for relativistic orbit')
plt.ylabel('y velocity (AU/yr)')
plt.xlabel('time (yr)')
plt.savefig('gr_yvel.png')

#plot angluar momentum over time
plt.figure(3)
plt.plot(t, L)
plt.grid()
plt.ylim([3.8,3.88])
plt.title('Angular momentum over time for relativistic orbit')
plt.xlabel('time (yr)')
plt.ylabel('angular momentum (M$_p$AU$^2$/yr)')
plt.savefig('gr_angular.png')

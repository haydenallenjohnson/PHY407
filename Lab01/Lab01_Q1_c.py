#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 11 21:27:01 2018

@author: Hayden
"""
'''
Code for part c) - numerical integration of equations of motion for
mercury under the influence of Newtonian gravity
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
    
    #calculate r^3 once first to save time
    r_cubed = (x[i]*x[i] + y[i]*y[i])**(1.5)
    
    #update velocities
    v_x[i+1] = v_x[i] - G*M_s*x[i]*dt/r_cubed
    v_y[i+1] = v_y[i] - G*M_s*y[i]*dt/r_cubed

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
plt.title('Newtonian orbit of Mercury')
plt.xlabel('x position (AU)')
plt.ylabel('y position (AU)')
plt.savefig('newton_orbit.png')

#plot x velocity vs. time
plt.figure(1)
plt.plot(t, v_x)
plt.grid()
plt.title('x velocity vs. time for Netownian orbit')
plt.ylabel('x velocity (AU/yr)')
plt.xlabel('time (yr)')
plt.savefig('newton_xvel.png')

#plot y velocity vs. time
plt.figure(2)
plt.plot(t, v_y)
plt.grid()
plt.title('y velocity vs. time for Netownian orbit')
plt.ylabel('y velocity (AU/yr)')
plt.xlabel('time (yr)')
plt.savefig('newton_yvel.png')

#plot angluar momentum over time
plt.figure(3)
plt.plot(t, L)
plt.grid()
plt.ylim([3.8,3.88])
plt.title('Angular momentum over time for Newtonian orbit')
plt.xlabel('time (yr)')
plt.ylabel('angular momentum (M$_p$AU$^2$/yr)')
plt.savefig('newton_angular.png')

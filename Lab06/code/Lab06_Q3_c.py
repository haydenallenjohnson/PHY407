#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 01:00:59 2018

@author: Hayden
"""

#import modules
import numpy as np
import matplotlib.pyplot as plt

#set values of constants
sigma = 1.0
epsilon = 1.0
m = 1.0

#specify number of interacting particles
N = 16

#create initial conditions (code adapted from handout)
Lx = 4.0
Ly = 4.0

dx = Lx/np.sqrt(N)
dy = Ly/np.sqrt(N)

x_grid = np.arange(dx/2, Lx, dx)
y_grid = np.arange(dy/2, Ly, dy)

xx_grid, yy_grid = np.meshgrid(x_grid, y_grid)

x_initial = xx_grid.flatten()
y_initial = yy_grid.flatten()

r_initial = np.zeros((N,2))
r_initial[:,0] = x_initial
r_initial[:,1] = y_initial

v_initial = np.zeros((N,2))

#define function to compute acceleration of N particles due to interactions
def accel(r):
    
    # r = [r1, r2, ... rN]
    # ri = [xi,yi]

    #create empty array to store accelerations
    a = np.zeros(r.shape)

    for i in range(len(r)):
        for j in range(i):
            
            #compute difference in positions
            x = r[i,0] - r[j,0]
            y = r[i,1] - r[j,1]
    
            #compute d^2
            d_sq = x*x + y*y
    
            #compute the value of the force
            force = 0.5*epsilon * (48*(sigma**12)*d_sq**(-7) - 24*(sigma**6)*d_sq**(-4))

            #compute the acceleration
            a_x = x*force/m
            a_y = y*force/m
            a_current = np.array([a_x,a_y])
            
            #add acceleration values to acceleration array
            a[i] += a_current
            a[j] += -a_current
            
    return a

#define function to compute acceleration of surround grid on center grid
def accel_grid(r_grid, r):
    
    # r_grid & r = [r1, r2, ... rN]; r_grid is one of the surrounding units
    # ri = [xi,yi]

    #create empty array to store accelerations
    a = np.zeros(r.shape)

    for i in range(len(r)):
        for j in range(len(r_grid)):
            
            #compute difference in positions
            x = r[i,0] - r_grid[j,0]
            y = r[i,1] - r_grid[j,1]
    
            #compute d^2
            d_sq = x*x + y*y
    
            #compute the value of the force
            if d_sq > 2:
                force = 0.5*epsilon * (48*(sigma**12)*d_sq**(-7) - 24*(sigma**6)*d_sq**(-4))
            else:
                force = 0
            #compute the acceleration
            a_x = x*force
            a_y = y*force
            a_current = np.array([a_x,a_y])
            
            #add acceleration values to acceleration array
            a[i] += a_current
            #a[j] += -a_current
            
    return a

#set time step and number of steps for solving ODEs
dt = 0.01
steps = 1000

#create array of time points
tpoints = np.arange(0,dt*(steps),dt)

#create 2d arrays to store position and velocity values for each particle
rpoints = np.zeros((N,2,steps))
vpoints = np.zeros((N,2,2*steps))

#initialize arays with initial values of position and velocity
rpoints[:,:,0] = r_initial
vpoints[:,:,0] = v_initial

#perform initialization of v(0.5h)
vpoints[:,:,1] = vpoints[:,:,0] + 0.5*dt*accel(rpoints[:,:,0])

#define a mod array
mod_array = np.full_like(np.zeros((N,2)), np.array([Lx,Ly]))

temp = []
temp_2 = []

#iterate over the number of steps specified
for i in range(steps-1):
    #execute equations 8-11 from handout to update position and velocity
    #adjust the position
    rpoints[:,:,i+1] = rpoints[:,:,i] + dt*vpoints[:,:,2*i+1]
    #mod the positions
    rpoints[:,:,i+1] = np.mod(rpoints[:,:,i+1], mod_array)
    #compute velocity change due to force within center grid
    k = dt*accel(rpoints[:,:,i+1])
    #add velocity change due to forces of surrounding grid
    """for j in range(-1,2,1):
        for m in range(-1,2,1):
            if m == 0 and j == 0:
                continue
            Lx = 4.0
            Ly = 4.0
            outside_x = rpoints[:,0,i+1]+Lx*j
            outside_y = rpoints[:,1,i+1]+Ly*m
            grid_r = np.array([outside_x,outside_y])
            temp.append(grid_r)
            outside_v = dt*accel_grid(grid_r, rpoints[:,:,i+1])
            temp_2.append(outside_v)
            k += outside_v
    """#update velocities
    vpoints[:,:,2*i+2] = vpoints[:,:,2*i+1] + 0.5*k
    vpoints[:,:,2*i+3] = vpoints[:,:,2*i+1] + k

plt.figure(0)
for i in range(N):
    plt.plot(rpoints[i,0,:],rpoints[i,1,:],'.')
    plt.savefig("../images/bdry.png", dpi=600)
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 17 10:35:31 2018

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
N = 2

#input initial positions and velocities
r1 = np.array([4.0,4.0])
r2 = np.array([5.2,4.0])
r_initial = np.array([r1,r2])

v1 = np.array([0.0,0.0])
v2 = np.array([0.0,0.0])
v_initial = np.array([v1,v2])

#define function to compute acceleration
def accel(r):
    
    # r = [[x1,y1],
    #      [x2,y2]]

    #compute difference in positions
    x = r[1,0] - r[0,0]
    y = r[1,1] - r[0,1]
    
    #compute d^2
    d_sq = x*x + y*y
    
    #compute the value of the force
    force = epsilon * (48*(sigma**12)*d_sq**(-7) - 24*(sigma**6)*d_sq**(-4))
    a_x = x*force/m
    a_y = y*force/m
    
    #create the arrays of acceletation
    a1 = np.array([a_x/2,a_y/2])
    a2 = -a1
    
    return np.array([-a1,-a2])

#set time step and number of steps for solving ODEs
dt = 0.01
steps = 100

#create array of time points
tpoints = np.arange(0,dt*(steps),0.5*dt)

#create 2d arrays to store position and velocity values for each particle
rpoints = np.zeros((N,2,steps))
vpoints = np.zeros((N,2,len(tpoints)))

#initialize arays with initial values of position and velocity
rpoints[:,:,0] = r_initial
vpoints[:,:,0] = v_initial

#perform initialization of v(0.5h)
vpoints[:,:,1] = vpoints[:,:,0] + 0.5*dt*accel(rpoints[:,:,0])

#iterate over the number of steps specified
for i in range(steps-1):
    #execute equations 8-11 from handout to update position and velocity
    rpoints[:,:,i+1] = rpoints[:,:,i] + dt*vpoints[:,:,2*i+1]
    k = dt*accel(rpoints[:,:,i+1])
    vpoints[:,:,2*i+2] = vpoints[:,:,2*i+1] + 0.5*k
    vpoints[:,:,2*i+3] = vpoints[:,:,2*i+1] + k


plt.figure(0)
#plt.plot(rpoints[0,0,:], rpoints[0,1,:], '.')
#plt.plot(rpoints[1,0,:], rpoints[1,1,:], '.')
plt.plot(rpoints[0,0,:],'.')
plt.plot(rpoints[1,0,:],'.')

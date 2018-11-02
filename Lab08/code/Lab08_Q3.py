# -*- coding: utf-8 -*-
"""
Solves Burger's equation using one Euler forward step and subsequent steps
using the leapfrog method. Plots the solution for various times.

Created on Thu Nov  1 07:36:40 2018
@author: Pierino Zindel
"""

#import libraries
import numpy as np
import pylab as plt


#define constants
epsilon = 1.0
dx = 0.02
dt = 0.005
Lx = 2*np.pi
Tf = 2

#compute space and time step numbers
Nx = int(Lx//dx)
Nt = int(Tf//dt)

#times for plotting
t1 = 0.0
t2 = 0.5
t3 = 1.0
t4 = 1.5

#offset for plot timing
delta = 1e-10

#create arrays
u = np.zeros([Nt, Nx]) #contains solutions throughout time
ut = np.empty(Nx, float) #contains solutions for a given time

#initial condition
#u(x, t=0) = sin(x)
x = np.linspace(0, Lx, Nx)
u[0] = np.sin(x)

#boundry conditions
#u(x=0, t) = 0
#u(x=Lx, t) = 0
ut[0] = 0
ut[-1] = 0
u[:,0] = 0
u[:,-1] = 0

#apply initial step with Euler forward
u[1] = np.sin(x)*(1 - epsilon*dt*np.cos(x))

#plot initial wave for t=0
plt.figure
plt.plot(x, u[0], label='t=0.0')

#loop over time, excluding the initial condition
for j in range(1,Nt-1):
    #current time
    t = j*dt
    #loop over x, excluding the boundry
    for i in range(1,Nx-1):
        #formula constant
        beta = epsilon*dt/dx
        #apply leapfrog step
        u[j+1,i] = u[j-1,i] - beta/2*((u[j,i+1])**2 - (u[j,i-1])**2)

    #plot the function at the given times
    if abs(t-t2) < delta:
        plt.plot(x, u[j], label='t=0.5')
    elif abs(t-t3) < delta:
        plt.plot(x, u[j], label='t=1.0')
    elif abs(t-t4) < delta:
        plt.plot(x, u[j], label='t=1.5')
        
#plotting values
plt.xlabel("Position, x")
plt.ylabel("Displacement, u(x,t)")
plt.title("Solutions to Burger's equation")
plt.grid()
plt.legend()
plt.savefig("../images/burgers.png", dpi=600)
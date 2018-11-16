# -*- coding: utf-8 -*-
"""
Solves Burger's equation using the Lax-Wendroff method. Plots the solutions at
times t=0, 0.5, 1.0, 1.5 seconds.

Created on Tue Nov 13 18:21:56 2018
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

#plot initial wave for t=0
plt.figure
plt.plot(x, u[0], label='t=0.0')

#loop over time and solve
for j in range(0,Nt-1):
    #current time
    t = j*dt
    #loop over x, excluding the boundry
    for i in range(1,Nx-1):
        #formula constant
        beta = epsilon*dt/dx
        #apply Lax-Wendoff method
        term_1 = u[j,i]
        term_2 = -(beta/4)*((u[j,i+1])**2 - (u[j,i-1])**2)
        term_3 = (beta**2/8)*(u[j,i] + u[j,i+1])*((u[j,i+1])**2 - (u[j,i])**2)
        term_4 = (beta**2/8)*(u[j,i] + u[j,i-1])*((u[j,i-1])**2 - (u[j,i])**2)
        u[j+1,i] = term_1 + term_2 + term_3 + term_4

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
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 14 10:53:07 2018

@author: Pierino
"""

import numpy as np
import pylab as plt
import dcst.py as dcst

#constant of pde
v = 100 #m/s
L = 1 #m
d = 0.1 #m
C = 1 #m/s
sigma = 0.3 #m
h = 1e-6 #s

#initial velocity profile
def psi0(x):
    return C*x*(L-x)*np.exp(-(x-d)**2/(2*sigma**2))/L**2

#grid points
N = 100
#grid spacing
a = L/N

#create arrays containing the initial conditions
x = np.linspace(0,L,N+1)
psi = psi0(x)
phi = np.zeros(N+1, float)
psi_step = np.copy(psi)
phi_step = np.copy(phi)

#loop time
t = 0.0
t_f = 0.101

#array of positions corresponding to each loop over time
traj = [phi]

#array of times corresponding to the list of trajectories
times = [t]

#frame count
frame = 0

#images save count
save = 0

#loop over time
while t < t_f:
    #constant of integration
    c = h*v**2/a**2
    #FTCS stepping
    phi_step[1:N] = phi[1:N] + h*psi[1:N]
    psi_step[1:N] = psi[1:N] + c*(phi[2:N+1] + phi[0:N-1] - 2*phi[1:N])
    
    #set old array to the new displacement
    phi = np.copy(phi_step)
    psi = np.copy(psi_step)    
    
    #increase time
    t += h
    
    #add current plot to container
    traj.append(np.copy(phi_step))
    times.append(t)
    
    #increase frame count for 'animation'
    frame += 1
    
    #plot frame
    if (frame % 500) == 0:
        frame = 0
        #plt.clf()
        plt.figure()
        plt.plot(x, phi)
        plt.xlim([0,1])
        plt.ylim([-0.0005, 0.0005])
        plt.xlabel("Position, $x$")
        plt.ylabel("Displacement, $\phi(x)$")
        plt.title("Piano string displacement at t=%.5f" %t)
        plt.grid()
        
        #save the plot for certain times
        if t > 0.1 and save == 2:
            plt.savefig("q2_t=100ms.png", dpi=600)
            save += 1
        elif t > 0.05 and save == 1:
            plt.savefig("q2_t=50ms.png", dpi=600)
            save += 1
        elif t > 0.005 and save == 0:
            plt.savefig("q2_t=2ms.png", dpi=600)
            save += 1
        
        plt.draw()
        plt.pause(0.01)



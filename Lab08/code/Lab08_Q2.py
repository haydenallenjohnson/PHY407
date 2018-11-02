# -*- coding: utf-8 -*-
"""
Created on Wed Oct 31 14:25:26 2018

@author: Pierino
"""

#import libraries
import numpy as np
import pylab as plt


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
N = 200
#grid spacing
a = L/N

x = np.linspace(0,L,N+1)
psi = psi0(x)
phi = np.zeros(N+1, float)
psi_step = psi
phi_step = phi

#loop
t = 0.0
t_f = 1.0

#array of positions corresponding to each loop over time
traj = [phi]

#image

#frame count
frame = 0

while t < t_f:
    #for i in range(1,N):
    c = h*v**2/a**2
    phi_step[1:N] = phi[1:N] + h*psi[1:N]
    psi_step[1:N] = psi[1:N] + c*(phi[2:N+1] + phi[0:N-1] - 2*phi[1:N])
    
    phi = phi_step
    psi = psi_step
    traj.append(phi_step)
    t += h
    
    frame += 1
    if (frame % 50) == 0:
        frame = 0
        plt.clf()
        plt.plot(x, phi)
        plt.xlim([0,1])
        plt.ylim([-0.0005, 0.0005])
        plt.draw()
        plt.pause(0.01)
    
#plt.plot(x, phi)
plt.xlabel("x")
plt.ylabel("phi")
plt.show()


# -*- coding: utf-8 -*-
"""
Solves the wave equation corresponding to the problem of a piano hammer
striking a piano string. Implements the spectral method and plots the 
displacement over time at times t = 2, 4, 6, 12, 100 milliseconds.

Created on Wed Nov 14 10:53:07 2018

@author: Pierino Zindel
"""

#import libraries
import numpy as np
import pylab as plt
import dcst as dcst

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
N = 101
#grid spacing
a = L/N

#create arrays containing the initial conditions
x = np.linspace(0,L,N)
psi_0 = psi0(x)
phi_0 = np.zeros(N, float)

#compute the coefficients of the intial state
psi_coeff = dcst.dst(psi_0)
phi_coeff = dcst.dst(phi_0)

#containers for the coefficients returned by the discrete sine transform
sol_coeff = []

#array of times to be plotted
t = np.array([2,4,6,12,100])/1000

#container for the solutions at the desired times
solutions = []

#loops over the desired times and compute the solution to the wave equation
for t_i in t:
    #compute the omega values for each coefficient
    k = np.arange(0,N,1)
    omega = np.pi*v*k/L
    #compute the coefficients of our solution
    sol_coeff = phi_coeff*np.cos(omega*t_i) + psi_coeff/omega*np.sin(omega*t_i)
    #compute the solution using the coefficients
    phi = dcst.idst(np.real(sol_coeff))
    solutions.append(phi)

#plotting fucntion
def plot_fig(x, phi, label, save, file_name):
    plt.figure()
    plt.plot(x, phi)
    plt.xlim([0,1])
    plt.ylim([-0.0005, 0.0005])
    plt.xlabel("Position, $x$")
    plt.ylabel("Displacement, $\phi(x)$")
    plt.title("Piano string displacement at " + label)
    plt.grid()
    if save:
        plt.savefig("../images/"+file_name, dpi=600)

#plot the solution for the desired times
plot_fig(x, solutions[0], 't=2 ms', True, 'wave_eqtn_2ms.png')
plot_fig(x, solutions[1], 't=4 ms', True, 'wave_eqtn_4ms.png')
plot_fig(x, solutions[2], 't=6 ms', True, 'wave_eqtn_6ms.png')
plot_fig(x, solutions[3], 't=12 ms', True, 'wave_eqtn_12ms.png')
plot_fig(x, solutions[4], 't=100 ms', True, 'wave_eqtn_100ms.png')


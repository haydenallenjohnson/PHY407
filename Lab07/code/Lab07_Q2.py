#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 24 10:03:24 2018

@author: Hayden
"""

#import modules
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as pc

def simpsons(x,f):
    '''
    computes the integral, using simpson's method, of a function given its
    values f evaluated on an array x
    '''
    N = len(x)
    h = (x[-1]-x[0])/N
    s = f[0] + f[-1]
    for k in range(1,N,2):
        s += 4*f[k]
    for k in range(2,N,2):
        s += 2*f[k]
    return h*s/3.0

# Constants
a = pc.physical_constants['Bohr radius'][0]
E0 = pc.physical_constants['Rydberg constant times hc in eV'][0] 
m = pc.m_e # electron mass
hbar = pc.hbar
e = pc.e # elementary charge
eps0 = pc.epsilon_0
pi = np.pi

#set n and l
n = 2
l = 0

#set step size and max r
h = 0.002*a
r_inf = 20*a

#define array of points r
r = np.arange(h,r_inf,h)

# Potential function
def V(x):
    return -e*e/(4*pi*eps0*x)

def f(val,x,E):
    R = val[0]
    S = val[1]
    dR = S/(x*x)
    dS = ((2*m*x*x/hbar**2)*(V(x)-E) + l*(l+1))*R
    return np.array([dR,dS],float)

# Calculate the wavefunction for a particular energy
def solve(E):
    R = 0.0
    S = 1.0
    vals = np.empty((len(r),2),float)
    vals[0] = np.array([R,S])
    
    for i in range(len(r)-1):
        k1 = h*f(vals[i],r[i],E)
        k2 = h*f(vals[i]+0.5*k1,r[i]+0.5*h,E)
        k3 = h*f(vals[i]+0.5*k2,r[i]+0.5*h,E)
        k4 = h*f(vals[i]+k3,r[i]+h,E)
        vals[i+1] = vals[i] + (k1+2*k2+2*k3+k4)/6

    return vals

# Main program to find the energy using the secant method
E1 = -15*e/n**2
E2 = -13*e/n**2
R2 = solve(E1)[-1,0]

target = e/1000
while abs(E1-E2) > target:
    R1,R2 = R2,solve(E2)[-1,0]
    E1,E2 = E2,E2-R2*(E2-E1)/(R2-R1)
    
psi = solve(E2)[:,0]

norm = simpsons(r,abs(psi)**2)

psi_normal = psi/np.sqrt(norm)
plt.figure(0)
plt.plot(r, psi_normal)
plt.grid()
plt.title('Normalized wave function for $n='+str(n)+', l='+str(l)+'$')
plt.ylabel('$\psi(r)$')
plt.xlabel('$r$ [m]')
plt.savefig('../images/q2_d_n='+str(n)+'_l='+str(l)+'.png')

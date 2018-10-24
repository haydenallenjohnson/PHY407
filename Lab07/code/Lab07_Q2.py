#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 24 10:03:24 2018

@author: Hayden
"""

#import modules
import numpy as np
import scipy.constants as pc

#define function for integration using simpson's rule
def simpsons(f,a,b,N):
    h = (b-a)/N
    s = f(a) + f(b)
    for k in range(1,N,2):
        s += 4*f(a+k*h)
    for k in range(2,N,2):
        s += 2*f(a+k*h)
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
n = 1
l = 0

#set step size and max r
h = 0.0002*a
r_inf = 20*a

# Potential function
def V(r):
    return -e*e/(4*pi*eps0*r)

def f(vals,r,E):
    R = vals[0]
    S = vals[1]
    dR = S/(r*r)
    dS = ((2*m*r*r/hbar**2)*(V(r)-E) + l*(l+1))*R
    return np.array([dR,dS],float)

# Calculate the wavefunction for a particular energy
def solve(E):
    R = 0.0
    S = 1.0
    vals = np.array([R,S],float)

    for r in np.arange(h,r_inf,h):
        k1 = h*f(vals,r,E)
        k2 = h*f(vals+0.5*k1,r+0.5*h,E)
        k3 = h*f(vals+0.5*k2,r+0.5*h,E)
        k4 = h*f(vals+k3,r+h,E)
        vals += (k1+2*k2+2*k3+k4)/6

    return vals[0]

# Main program to find the energy using the secant method
E1 = -15*e/n**2
E2 = -13*e/n**2
R2 = solve(E1)

target = e/1000
while abs(E1-E2)>target:
    R1,R2 = R2,solve(E2)
    E1,E2 = E2,E2-R2*(E2-E1)/(R2-R1)

print("E =",E2/e,"eV")
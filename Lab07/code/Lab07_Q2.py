#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 24 10:03:24 2018

@author: Hayden

Program which uses the shooting and RK4 methods to find the energy eigenvalue
and corresponding eigenfunction for the radial part of the Schrodinger
equation with the specified values of n and l, and plots these eigenfunctions 
along with their analytic counterparts on the same plot
"""

#import modules
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as pc

#define constants
a = pc.physical_constants['Bohr radius'][0]
E0 = pc.physical_constants['Rydberg constant times hc in eV'][0] 
m = pc.m_e #electron mass
hbar = pc.hbar
e = pc.e #elementary charge
eps0 = pc.epsilon_0
pi = np.pi

#define integration function
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

#define function to compute the analytic solutions R(r):
def R_theory(x, n, l):
    if n==1:
        return 2*np.exp(-x/a)/(a**1.5)
    
    elif n==2 and l==0:
        return (1/(2*np.sqrt(2)*a**1.5))*(2-(x/a))*np.exp(-x/(2*a))
    
    elif n==2 and l==1:
        return (1/(2*np.sqrt(6)*a**1.5))*(x/a)*np.exp(-x/(2*a))
    
#set n and l
n = 2
l = 1

#set step size and max r
h = 0.0002*a
r_inf = 20*a

#define array of points r
r = np.arange(h,r_inf,h)

#define potential function
def V(x):
    return -e*e/(4*pi*eps0*x)

#define function to output values of derivatives of coupled ODE system
def f(val,x,E):
    R = val[0]
    S = val[1]
    dR = S/(x*x)
    dS = ((2*m*x*x/hbar**2)*(V(x)-E) + l*(l+1))*R
    return np.array([dR,dS],float)

#calculate the wavefunction for a particular energy
def wave(E):
    R = 0.0
    S = 1.0
    vals = np.empty((len(r),2),float)
    vals[0] = np.array([R,S])
    
    for i in range(len(r)-1):
        k1 = h*f(vals[i], r[i], E)
        k2 = h*f(vals[i]+0.5*k1, r[i]+0.5*h, E)
        k3 = h*f(vals[i]+0.5*k2, r[i]+0.5*h, E)
        k4 = h*f(vals[i]+k3, r[i]+h, E)
        vals[i+1] = vals[i] + (k1+2*k2+2*k3+k4)/6

    return vals

# Main program to find the energy using the secant method
E1 = -15*e/n**2
E2 = -13*e/n**2
R2 = wave(E1)[-1,0]

#specify target accuracy
target = e/1000

#iteratively apply the secant method until the target accuracy is reached
while abs(E1-E2) > target:
    R1, R2 = R2, wave(E2)[-1,0]
    E1, E2 = E2, E2-R2*(E2-E1)/(R2-R1)
    
print('Calculated: ',E2/e)
print('Theory: ',-E0/n**2)

#compute the wave function for the calculated energy
psi = wave(E2)[:,0]

#compute the norm of psi and normalize the wavefunction
norm = simpsons(r,abs(psi)**2)
psi_normal = psi/np.sqrt(norm)

#compute the theoretical wavefunction
psi_theory = R_theory(r, n, l)

#normalize the theoretical wave function
norm_theory = simpsons(r,abs(psi_theory)**2)
psi_theory_normal = psi_theory/np.sqrt(norm_theory)

#plot normalized wavefunctions (calculated and theoretical)
plt.figure(0)
plt.plot(r, psi_normal, label='Computed')
plt.plot(r, psi_theory_normal, '--', label='Theoretical')
plt.grid()
plt.title('Normalized wave function for $n='+str(n)+', l='+str(l)+'$')
plt.ylabel('$R(r)$')
plt.xlabel('$r$ [m]')
plt.legend(loc='upper right')
plt.savefig('../images/q2_d_n='+str(n)+'_l='+str(l)+'.png')

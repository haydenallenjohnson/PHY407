# -*- coding: utf-8 -*-
"""
Computes the eigenvalues, eigenvectors, wavefunction and probability density
for the asymmetric quantum well.

Created on Fri Oct  5 12:57:26 2018
@author: Pierino Zindel
"""

import scipy.constants as const
import numpy as np

#define constants
L = 5e-10 #m
a = 10 #eV
M = const.electron_mass #kg
e = const.elementary_charge #coulomb
pi = const.pi
hbar = const.hbar
eV = const.eV #joules

def matrix(N):
    matrix = np.zeros((N,N))
    for i in range(N):
        m = i+1
        for j in range(N):
            n = j+1
            if m == n:
                matrix[i,j] = 0.5*a + pi*pi*hbar*hbar*m*m/(2*M*L*L)/eV
            elif (m % 2 == 0 and n % 2 == 0) or (m % 2 == 1 and n % 2 == 1):
                matrix[i,j] = 0
            else:
                matrix[i,j] = -8*a*m*n/(pi*pi*((m*m-n*n)**2))
    return matrix

H_1 = matrix(N=10)
evals_1, evect_1 = np.linalg.eigh(H_1)

H_2 = matrix(N=100)
evals_2, evect_2 = np.linalg.eigh(H_2)

print(evals_1, '\n')
print(evals_2[:10], '\n')



#Define wavefunction function
def wavefunc(evect,x):
    value = 0
    for n in range(len(evect)):
        value += evect[n]*np.sin(pi*(n+1)*x/L)
    return value

#define range of x from 0 to L
x = np.linspace(0,L,100)

#compute wavefunctions
wave_1 = wavefunc(evect_1[0], x) 
wave_2 = wavefunc(evect_1[1], x)
wave_3 = wavefunc(evect_1[2], x)





# =============================================================================
# 
# #define integration function
# def integ(y):
#     return abs(y)**2
# 
# #Number of integration slices
# N = 1000
# #Integraing interval
# a = 0.0
# b = L
# #Width of the intergrating bins
# h = ((b) - (a))/N
# ##Sum the integrand using the trapazoid integration method
# I = h*(0.5*integ(a) + 0.5*integ(b))
# for k in range(1,N):
#     I += h*(integ(a + k*h))
# =============================================================================


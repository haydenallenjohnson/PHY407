# -*- coding: utf-8 -*-
"""
Computes the eigenvalues, eigenvectors, wavefunction and probability density
for the asymmetric quantum well.

Created on Fri Oct  5 12:57:26 2018
@author: Pierino Zindel
"""

#import modules
import scipy.constants as const
import numpy as np
import pylab as plt

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
evals = np.linalg.eigvalsh(H_1)
evals = np

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

#define integration function
def integ(x,evect):
    return abs(wavefunc(evect,x))**2

#Number of integration slices
N = 1000

#Integraing interval
a = 0.0
b = L

#Width of the intergrating bins
h = (b-a)/N

#Sum the probability density of ground state using the trapazoid integration method
I_1 = h*(0.5*integ(a, evect_1[0]) + 0.5*integ(b, evect_1[0]))
for k in range(1,N):
    I_1 += h*(integ(a + k*h, evect_1[0]))
    
##Sum the probability density of 1st excited state using the trapazoid integration method
I_2 = h*(0.5*integ(a, evect_1[1]) + 0.5*integ(b, evect_1[1]))
for k in range(1,N):
    I_2 += h*(integ(a + k*h, evect_1[1]))
    
##Sum the probability density of 2nd excited state using the trapazoid integration method
I_3 = h*(0.5*integ(a, evect_1[2]) + 0.5*integ(b, evect_1[2]))
for k in range(1,N):
    I_3 += h*(integ(a + k*h, evect_1[2]))

#divide values of the 
wave_1 = wave_1/I_1
wave_2 = wave_2/I_2
wave_3 = wave_3/I_3

plt.figure()
plt.grid()
plt.ylabel("$|\psi(x)|^2$")
plt.xlabel("x")
plt.suptitle("Asymmetric Quantum Well Probability Densities")
plt.plot(x,abs(wave_1)**2,label="ground state")
plt.plot(x,abs(wave_2)**2,label="1st excited state")
plt.plot(x,abs(wave_3)**2,label="2nd excited state")
plt.legend()

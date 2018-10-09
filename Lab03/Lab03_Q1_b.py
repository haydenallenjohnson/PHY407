#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 26 10:38:01 2018

@author: Hayden
"""

#import functions
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf
from Lab03_Q1_functions import trapezoid, simpsons, gauss, f
from gaussxw import gaussxwab

#define endpoints
a = 0
b = 3

#create array with values of N1
N = np.array([10,30,100,300,1000,3000,10000])

#create arrays to store values of integral computations
trap_I1 = np.zeros(len(N))
trap_I2 = np.zeros(len(N))
simp_I1 = np.zeros(len(N))
simp_I2 = np.zeros(len(N))
gauss_I1 = np.zeros(len(N))
gauss_I2 = np.zeros(len(N))

#iterate over values of N and compute I1 and I2 for all methods
for i in range(len(N)):
    trap_I1[i] = trapezoid(f,a,b,N[i])
    trap_I2[i] = trapezoid(f,a,b,2*N[i])
    simp_I1[i] = simpsons(f,a,b,N[i])
    simp_I2[i] = simpsons(f,a,b,2*N[i])

    x1,w1 = gaussxwab(N[i], a, b)
    x2,w2 = gaussxwab(2*N[i], a, b)
    gauss_I1[i] = gauss(f,x1,w1)
    gauss_I2[i] = gauss(f,x2,w2)

#calculate error estimates
trap_eps = abs(trap_I2 - trap_I1)/3.0
simp_eps = abs(simp_I2 - simp_I1)/15.0
gauss_eps = abs(gauss_I2 - gauss_I1)

#calculate actual errors
trap_err = abs(trap_I2 - erf(3))
simp_err = abs(simp_I2 - erf(3))
gauss_err = abs(gauss_I1 - erf(3))

#plot estimates and actual errors
plt.figure(0)
plt.xscale('log')
plt.yscale('log')
plt.plot(N, trap_eps, label='Trapezoid estimate')
plt.plot(N, simp_eps, label='Simpson\'s estimate')
plt.plot(N, gauss_eps, label='Gaussian estimate')
plt.plot(N, trap_err, '--', label='Trapezoid error')
plt.plot(N, simp_err, '--', label='Simpson\'s error')
plt.plot(N, gauss_err, '--', label='Gaussian error')
plt.legend(loc='upper right', fontsize=8)
plt.grid()
plt.title('Errors calculating erf(3) with different integration methods')
plt.xlabel('N')
plt.ylabel('Error')
plt.savefig('graph_1b.png')
plt.show()

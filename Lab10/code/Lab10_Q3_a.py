#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 17:27:22 2018

@author: Hayden

Program which computes the integral of the specified integrand function using 
the Monte Carlo mean value and importance sampling methods. For each method the
integral is computed n=100 times and the computed values are plotted in a
histogram.
"""

#import modules
import numpy as np
import matplotlib.pyplot as plt

#define integrand function
def f(x):
    return 1/(np.sqrt(x)*(1+np.exp(x)))

#define weighting function
def w(x):
    return 1/np.sqrt(x)

#define parameters
N = 10000 #number of sample points
n = 100 #number of times to repeat the computation

#specify bounds of integration
a = 0
b = 1

#create empty arrays to store values of sums
s_mv = np.empty(n)
s_is = np.empty(n)

#iterate over the number of computations
for i in range(n):
    #create uniform array of x points to sample at for mean value method
    x_mv = np.random.rand(N)*(b-a) + a
    
    #calculate sum for mean value method
    s_mv[i] = np.sum(f(x_mv))
    
    #create non-uniform array of x points to sample at for importance sampling
    x_is = (np.random.rand(N)**2)*(b-a) + a
    
    #compute sum for importance sampling method
    s_is[i] = np.sum(f(x_is)/w(x_is))
    
#compute integrals from the sums
I_mv = (b-a)*s_mv/N
I_is = 2*s_is/N

#create histograms
plt.figure(0)
plt.hist(I_mv, 10, range=[0.8, 0.88])
plt.title('Mean value method')
plt.xlabel('result')
plt.ylabel('frequency (counts)')
plt.savefig('../images/q3_a_mv.png')

plt.figure(1)
plt.hist(I_is, 10, range=[0.8, 0.88])
plt.title('Importance sampling method')
plt.xlabel('result')
plt.ylabel('frequency (counts)')
plt.savefig('../images/q3_a_is.png')


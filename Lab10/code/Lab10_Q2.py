#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 14:03:20 2018

@author: Hayden
"""

#import modules
import numpy as np

#define parameters
N = 1000000 #number of points
n = 10 #number of dimensions

#compute volume of rectangle of integration
vol = 2**n

#define indicator function f to be integrated
def f(x):
    '''
    f takes value 1 when inside the sphere and 0 when outside
    '''
    r = np.sqrt(np.sum(x*x))
    if r <= 1:
        return 1
    elif r > 1:
        return 0

#define array of N x values where each component is a random number generated
#from a uniform distribution over (-1,1)
x_array = np.random.rand(N,n)*2 - 1

#calculate sum of f(x), and f(x)^2, for each point x in the array
s_f = 0
s_f_sq = 0
for i in range(N):
    val = f(x_array[i])
    s_f += val
    s_f_sq += val*val

#compute average vaule of f and f^2
f_av = s_f/N
f_sq_av = s_f_sq/N

#compute the value of the integral
I = vol*f_av

#compute error in integral estimate
var_f = f_sq_av - f_av*f_av
error = vol*np.sqrt(var_f)/np.sqrt(N)

#print results
print('Value of integral:',I)
print('Error:',error)
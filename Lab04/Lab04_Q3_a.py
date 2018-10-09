#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  3 09:55:39 2018

@author: Hayden
"""

#import modules
import numpy as np
import matplotlib.pyplot as plt

c = np.arange(0,3.01,0.01)

x_final = np.zeros(len(c))

def f(x, c):
    return 1 - np.exp(-c*x)

def f_1(x, c):
    return c*np.exp(-c*x)

def error(x1,x2,c):
    num = x1-x2
    a = f_1(x1,c)
    if a==0:
        return 1e-8
    elif a==1:
        return 1
    else:
        denom = 1 - (1/a)
        return abs(num/denom)

for i in range(len(c)):
    
    x_old = 1
    x_new = f(x_old, c[i])

    while error(x_old, x_new, c[i]) > 1e-6:
        x_old = x_new
        x_new = f(x_new, c[i])
    
    x_final[i] = x_new
    
plt.figure(0)
plt.plot(c, x_final)
plt.xlim([0,3])
plt.xlabel('c')
plt.ylabel('x')
plt.title('Solution to $x = 1-e^{-cx}$ as a function of c')
plt.savefig('3a_graph.png')
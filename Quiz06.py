#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 14 12:40:58 2018

@author: Hayden
"""

import numpy as np

#define function to use 2nd order runge-kutta method
def rk(f, x0, t0, t1, N):
    h = (t1-t0)/N
    
    x = x0
    t = t0
    
    for i in range(N):
        k1 = h*f(x,t)
        k2 = h*f(x + 0.5*k1, t + 0.5*h)
        x += k2
        t += h
        
    return x, t

def f(x, t):
    return -x**4 + np.cos(t)

x_final, t = rk(f, 0, 0, 2, 10)

print(x_final)
print (t)
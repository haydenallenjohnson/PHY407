#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 25 22:47:09 2018

@author: Hayden
"""

#import modules
import numpy as np

#define function for integration using trapezoid rule
def trapezoid(f,a,b,N):
    h = (b-a)/N
    s = 0.5*f(a) + 0.5*f(b)
    for k in range(1,N):
        s += f(a+k*h)
    return h*s

#define function for integration using simpson's rule
def simpsons(f,a,b,N):
    h = (b-a)/N
    s = f(a) + f(b)
    for k in range(1,N,2):
        s += 4*f(a+k*h)
    for k in range(2,N,2):
        s += 2*f(a+k*h)
    return h*s/3.0

#define function for integration using gaussian quadrature
def gauss(f,x,w):
    s = 0.0
    for i in range(len(x)):
        s += f(x[i])*w[i]
    return s

#define integrand function for part ai
def f(x):
    return (2.0/np.sqrt(np.pi))*np.exp(-x*x)

#define derivative of f
def f_1(x):
    return -2*x*f(x)

#define third deriviative of f
def f_3(x):
    return (3 - 2*x*x)*4*x*f(x)

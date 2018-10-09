#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 20 14:45:18 2018

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

#define integrand function for part ai
def f(x):
    return (2.0/np.sqrt(np.pi))*np.exp(-x*x)

#define derivative of f
def f_1(x):
    return -2*x*f(x)

#define third deriviative of f
def f_3(x):
    return (3 - 2*x*x)*4*x*f(x)

#define function to compute bessel function
def J(m,x):
    #define the integrand function f
    def f(theta):
        return np.cos(m*theta - x*np.sin(theta))

    #set up parameters for simpson's method integration
    a = 0
    b = np.pi
    N = 1000
    
    #compute value of integral
    integral = simpsons(f, a, b, N)
    
    #return result
    return integral/np.pi

#define function to compute intensity as a function of radius and wavelength
def I(r, lam):
    k = 2*np.pi/lam   
    frac = J(1,k*r)/(k*r)
    return frac*frac

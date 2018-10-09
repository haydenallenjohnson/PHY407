#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 25 22:59:08 2018

@author: Hayden
"""

#import functions
from scipy.special import erf
from Lab03_Q1_functions import trapezoid, simpsons, gauss, f
from gaussxw import gaussxwab

#set parameters of integration
N = 1000
a = 0
b = 3

#compute output 
erf_out = erf(3)
trap_out = trapezoid(f, a, b, N)
simp_out = simpsons(f, a, b, N)

x,w = gaussxwab(N, a, b)
gauss_out = gauss(f, x, w)

print('scipy.erf: ' + str(erf_out))
print('Trapezoid rule: ' + str(trap_out))
print('Spimson\'s rule: ' + str(simp_out))
print('Guassian quad: ' + str(gauss_out))

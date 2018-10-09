#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 19 11:47:59 2018

@author: Hayden
"""

#import modules
from scipy.special import erf
from Lab02_Q2_functions import trapezoid, simpsons, f, f_1, f_3

#create parameters for integration
a = 0
b = 3
N1 = 10
N2 = 2*N1
h = (b-a)/N2

#compute scipy.special.erf
erf_2 = erf(3)

#compute integrals with trapezoid rule for both Ns
trap_1 = trapezoid(f, a, b, N1)
trap_2 = trapezoid(f, a, b, N2)

#compute integrals with simpson's rule for both Ns
simp_1 = simpsons(f, a, b, N1)
simp_2 = simpsons(f, a, b, N2)

#compute estimated error for trapezoid method
trap_eps2 = abs(trap_2 - trap_1)/3.0

#compute esimated error for simpson's method
simp_eps2 = abs(simp_2 - simp_1)/15.0

#compute errors from euler-maclaurin formula
trap_em_eps = h*h*(f_1(a) - f_1(b))/12.0
simp_em_eps = h*h*h*h*(f_3(a) - f_3(b))/180.0

#print estimated errors
print('Trapezoid estimated error: ' + str(trap_eps2))
print('Trapezoid Eulier-Maclaurin error: ' + str(trap_em_eps))
print('Simposon\'s estimated error: ' + str(simp_eps2))
print('Simposon\'s Eulier-Maclaurin error: ' + str(simp_em_eps))
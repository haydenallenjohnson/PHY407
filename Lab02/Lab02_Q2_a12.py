#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 19 09:13:29 2018

@author: Hayden
"""

#import modules
from scipy.special import erf
from time import time
from Lab02_Q2_functions import trapezoid, simpsons, f 

#set parameters of integration
N = 10000
a = 0
b = 3

#compute output of scipy.special.erf(3)
t1 = time()
erf_out = erf(3)
t2 = time()
t_erf = t2 - t1

#compute integral using trapezoid method
t1 = time()
trap_out = trapezoid(f, a, b, N)
t2 = time()
t_trap = t2 - t1

#compute integral using simpson's method
t1= time()
simp_out = simpsons(f, a, b, N)
t2 = time()
t_simp = t2 - t1

#calculate relative errors for trapezoid and simpson's rule
trap_error = abs(erf_out - trap_out)/erf_out
simp_error = abs(erf_out - simp_out)/erf_out


#print output
#print('scipy.erf: ' + str(erf_out))
print('scipy.erf time: ' + str(t_erf))
#print('Trapezoid rule: ' + str(trap_out))
print('Trapezoid error: ' + str(trap_error))
print('Trapezoid time: ' + str(t_trap))
#print('Spimson\'s rule: ' + str(simp_out))
print('Simpson\'s error: ' + str(simp_error))
print('Simpson\'s time: ' + str(t_simp))

# -*- coding: utf-8 -*-
"""
Computes the integral of the formula for the total energy per unit area 
radiated by a black body, and computes the value of the Stefan-Boltzmann 
constant.

Created on Wed Sep 19 11:00:28 2018
@author: Pierino Zindel
"""
#Import libraries
import numpy as np
import scipy.constants as const

#Define function to be integrated over
def func(z):
    return (z**3/(((1-z)**5)*(np.exp(z/(1-z))-1)))

#Number of integration slices
N = 1000
#Interval offset to avoid division by zero within the integration
delta = 1e-10
#Integraing interval
a = 0.0
b = 1.0
#Width of the intergrating bins
h = ((b-delta) - (a+delta))/N

#Sum the integrand using the trapazoid integration method
I = h*(0.5*func(a+delta) + 0.5*func(b-delta))
for k in range(1,N):
    I += h*(func(a+delta + k*h))

print("The integral is equal to: ", I)

#Compute Stefan-Boltzmann constant using our computed integrand
sb_const = (const.k**4 * I)/(4 * np.pi**2 * const.c**2 * const.hbar**3)

#Compare Stefan-Boltzmann constant defined in scipy
print("The Stefan-Boltzmann constant computed using the integral is: ", 
      sb_const)
print("The Stefan-Boltzmann constant as defined in the scipy module is: ", 
      const.Stefan_Boltzmann)


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 19 23:28:18 2018

@author: Hayden
"""

#import modules
import numpy as np
import matplotlib.pyplot as plt
from Lab02_Q2_functions import J
from scipy.special import jv

#create x array
x = np.arange(0, 20, 0.1)

#create arrays of function values
J0 = J(0,x)
J1 = J(1,x)
J2 = J(2,x)

J0_sp = jv(0,x)
J1_sp = jv(1,x)
J2_sp = jv(2,x)

#plot functions
plt.figure(0)
plt.plot(x, J0, label='simpson\'s $J_0$')
plt.plot(x, J0_sp, linestyle='--', label='scipy $J_0$')
plt.plot(x, J1, label='simpson\'s $J_1$')
plt.plot(x, J1_sp, linestyle='--', label='scipy $J_1$')
plt.plot(x, J2, label='simpson\'s $J_2$')
plt.plot(x, J2_sp, linestyle='--', label='scipy $J_2$')
plt.grid()
plt.xlim([0,20])
plt.legend(loc='upper right')
plt.xlabel('$x$')
plt.ylabel('$J(x)$')
plt.title('Plot of Bessel functions $J_0$, $J_1$ and $J_2$')
plt.savefig('bessel_graph_compare.png')
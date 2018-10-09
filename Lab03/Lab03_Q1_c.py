#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 26 14:59:01 2018

@author: Hayden

Code computes probability of blowing snow as a function of the average
air temperature, wind speed, and surface snow age, and plots probability
as a function of temperature for different values of the other
parameters
"""

#import modules
import numpy as np
import matplotlib.pyplot as plt
from Lab03_Q1_functions import gauss, f
from gaussxw import gaussxw

#define function to compute U_bar
def u_bar(T_a, t_h):
    return 11.2 + 0.365*T_a + 0.00706*T_a*T_a +0.9*np.log(t_h)

#define function to compute delta
def delta(T_a):
    return 4.3 + 0.145*T_a + 0.00196*T_a*T_a

#define integrand function g which is 0.5*f
def g(x):
    return 0.5*f(x)

#define parameter N
N = 100

#calculate points and weights for gaussian integration with N points
x,w = gaussxw(N)

#define function to compute probability of blowing snow
def P(u_10, T_a, t_h):
    
    #compute values of u_bar and delta
    u_bar_var = u_bar(T_a, t_h)
    delta_var = delta(T_a)
    
    #compute bounds of integration
    t_0 = (u_bar_var - u_10)/(np.sqrt(2)*delta_var)
    t_1 = u_bar_var/(np.sqrt(2)*delta_var)
    
    #apply rescaling to x and w based on endpoints
    x_current = 0.5*(t_1-t_0)*x+0.5*(t_1+t_0)
    w_current = 0.5*(t_1-t_0)*w
    
    return gauss(g,x_current,w_current)
 
#vectorize P function
P_vec = np.vectorize(P)
    
#define values of u_10 and t_h
u_10_vals = (6,8,10)
t_h_vals = (24,48,72)

#create tuples of colors and linestyles
colors = ('b','g','r')
lines = ('-','--',':')

#create array of values for T_a
T_a_array = np.arange(-50,10,0.1)

#create plot
plt.figure(0)

#iterate over values of u_10 and t_h
for (u_10_val, color) in zip(u_10_vals, colors):
    for (t_h_val, line) in zip(t_h_vals, lines):
        
        #create strings
        plot_str = color + line
        label_str = '$u_{10}$='+str(u_10_val)+', $t_h$='+str(t_h_val)
        
        #calculate P values and plot them
        P_array = P_vec(u_10_val, T_a_array, t_h_val)        
        plt.plot(T_a_array, P_array, plot_str, label=label_str)
        
#format plot
plt.title('Probability of blowing snow as a function of temperature')
plt.ylabel('Probability of blowing snow')
plt.xlabel('Temperature ($^\circ$C)')
plt.legend(loc='upper right', fontsize=7)
plt.xlim([-50,10])
plt.grid()
plt.savefig('graph_1c.png')
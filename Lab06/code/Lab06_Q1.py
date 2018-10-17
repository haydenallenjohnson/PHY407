# -*- coding: utf-8 -*-
"""
Numerically solve two second order ode corresponding to a ball bearing orbiting
aroung a rod in space. Decomposes the ode to 4 first order odes and solves 
using 4th-order Runge-Kutta method.

Created on Sun Oct 14 10:14:04 2018
@author: Pierino Zindel
"""

#import libraries
import numpy as np
import pylab as plt

#define ode function of 4 first order odes
def func(r, t):
    #unpack initial conditions
    x, y, vx, vy = r
    
    #constants
    G = 1.0
    M = 10.0
    L = 2.0
    
    #compute new conditions
    R = (x**2 + y**2)
    coeff = -G*M/R/(R + L**2/4)**0.5
    
    Dvx = x*coeff
    Dvy = y*coeff
    Dx = vx
    Dy = vy
    
    return np.array([Dx,Dy,Dvx,Dvy])

#define RK4 method function
def RK4(h, f, r, t):
    k1 = h*f(r,t)
    k2 = h*f(r + k1/2, t + h/2)
    k3 = h*f(r + k2/2, t + h/2)
    k4 = h*f(r + k3, t + h)    
    return r + (k1 + 2*k2 + 2*k3 + k4)/6

#ode solve parameters
a = 0.0
b = 10.0
N = 1000
h = (b-a)/N

#initial conditions & containers
r = np.array([1.0, 0.0, 0.0, 1.0]) #[x,y,vx,vy]
solution = []
tpoints = np.arange(a,b,h)

#solve ode
for t in tpoints:
    solution.append(r)
    r = RK4(h, func, r, t)

#extract x&y values from solution container
solution = np.array(solution)
x = solution[:,0]
y = solution[:,1]

#plot resulting function
plt.plot(x, y)
plt.xlabel("x")
plt.ylabel("y")
plt.title("Orbit of ball bearing around rod in space")
plt.grid()
plt.savefig("../images/q1_orbit.png", dpi=600)


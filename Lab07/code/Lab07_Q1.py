# -*- coding: utf-8 -*-
"""
Numerically solve two second order ode corresponding to a ball bearing orbiting
aroung a rod in space. Decomposes the ode to 4 first order odes and solves 
using 4th-order Runge-Kutta method with adaptive time step. Compares plots 
and computation time of using adaptive time step and using fixed time step.

Created on Wed Oct 24 10:49:54 2018
@author: Pierino Zindel
"""

#import libraries
import numpy as np
import pylab as plt
from time import time

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

#calculate accuracy value rho
def rho_calc(h, delta, x1, x2, y1, y2):
    #compute the Euclidean error
    eps_x = (x1 - x2)/30
    eps_y = (y1 - y2)/30
    div = np.sqrt((eps_x)**2 + (eps_y)**2)
    
    #raise error if divisor is invalid
    if div == 0:
        raise Exception("Euclidean error is 0; division by zero in rho_cal")
    return h*delta/div

########################################################################
#code copied from Lab06 Q1 to plot trajectories using constant time step

#define number of steps and step size
N = 1000
h = 0.01

#initial conditions
r = np.array([1.0, 0.0, 0.0, 1.0]) #[x,y,vx,vy]
#container for solution values
sol = []
#array of times
tpoints = np.arange(0,N*h,h)

#record start time
start = time()

#solve ode
for t in tpoints:
    sol.append(r)
    r = RK4(h, func, r, t)

#record stop time
stop = time()

#print the computation time
print("Elapsed time to solve with constant time step: ", stop-start)
print("Number of steps N =", N, " of size h =", h)

#extract x&y values from solution container
sol = np.array(sol)
x = sol[:,0]
y = sol[:,1]

#plot resulting trajectory
plt.figure()
plt.plot(x, y, label='constant step solution')
plt.grid()

########################################################################

#error per time step
delta = 1e-6
#number of steps
N = 1000
#initial step size
h = 0.01
#initial conditions
r0 = np.array([1.0, 0.0, 0.0, 1.0]) #[x,y,vx,vy]
#define container for solution
sol = [r0]
#keep track of times for stepping
time_1 = 0
time_2 = time_1 + h
#initialize array of steps
steps = np.arange(N//2)
#initialize containers for the steps sizes and the times they are taken
time_step = []
timing = []
#record start time
start = time()

#solve the ode
for i in steps:
    time_step.append(h)
    time_step.append(2*h)
    timing.append(time_1)
    timing.append(time_1+h)
    #initialize rho
    rho = 0
    
    #compute and adjust time steps until accurate enough
    while rho < 1:
        #take two single steps
        r_one_step = RK4(h, func, r0, time_1)
        r_two_step = RK4(h, func, r_one_step, time_2)
        #take one double step
        r_double_step = RK4(2*h, func, r0, time_1)
        #record the estimated positions
        x1 = r_two_step[0]
        x2 = r_double_step[0]
        y1 = r_two_step[1]
        y2 = r_double_step[1]
        #compute rho and new time step
        rho = rho_calc(h, delta, x1, x2, y1, y2)
        h = h*rho**0.25
        
    #store the steps
    sol.append(r_one_step)
    sol.append(r_two_step)
    #shift the time
    time_1 += 2*h
    time_2 = time_1 + h
    #shift starting r values
    r0 = r_two_step
    
#record stop time
stop = time()

#print the computation time
print("Elapsed time to solve with adaptive time step: ", stop-start)
print("Number of steps N =", N, "of initial size h = 0.01")

#extract x&y values from solution container
sol = np.array(sol)
x = sol[:,0]
y = sol[:,1]

#plot resulting function
plt.plot(x, y, 'k.', label='adaptive time steps')
plt.xlabel("x")
plt.ylabel("y")
plt.title("Orbit of ball bearing around rod in space")
plt.grid()
plt.savefig("../images/q1_orbit.png", dpi=600)

#plot the step size as a function of time
plt.figure()
plt.plot(timing, time_step, '-')
plt.xlabel("Time, t (s)")
plt.ylabel("Time step, h (s)")
plt.title("Time steps over time corresponding to orbit")
plt.grid()
plt.savefig("../images/q1_orbit_step.png", dpi=600)


# -*- coding: utf-8 -*-
"""
Created on Tue Oct 16 12:54:17 2018
@author: Pierino
"""

#import libraries
import numpy as np
import pylab as plt

#function of equations of motion
def func(r, v, t):
    #unpack initial conditions
    x, y = r
    vx, vy = v
    
    #constants
    epsilon = 1.0
    m = 1.0
    sigma = 1.0
    
    #compute new conditions
    R = (x**2 + y**2)
    c1 = (64.0 * epsilon * R**(-9) * sigma**(12))/m
    c2 = (-24.0* epsilon * R**(-4) * sigma**(6))/m
    
    Dvx = (c1+c2)*x/2
    Dvy = (c1+c2)*y/2
    Dx = vx
    Dy = vy
    
    return np.array([Dx,Dy,Dvx,Dvy])

#define RK4 method function
#def RK4(h, f, r, t):
#    k1 = h*f(r,t)
#    k2 = h*f(r + k1/2, t + h/2)
#    k3 = h*f(r + k2/2, t + h/2)
#    k4 = h*f(r + k3, t + h)    
#    return r + (k1 + 2*k2 + 2*k3 + k4)/6


#def verlet():
#    v(t+h/2)=v(t) +h*f(r(t),t)/2
#    
#    #repeat
#    r(t+h)=r(t) +h*v(t+h/2)
#    k = h*f(r(t+h), t+h)
#    v(t+h)=v(t+h/2)+k/2
#    v(t+3*h/2)=v(t+h/2)+k
    
def solve(r1, r2, f, h, tpoints):
    N = int(len(tpoints)/2)
    s1 = np.zeros([N, 2])
    s2 = np.zeros([N, 2])
    v1 = np.zeros([N, 2])
    v2 = np.zeros([N, 2])
    
    s1[0][0] = r1[0]
    s1[0][1] = r1[1]
    s2[0][0] = r2[0]
    s2[0][1] = r2[1]
    v1[0][0] = r1[2]
    v1[0][1] = r2[3]
    v2[0][0] = r2[2]
    v2[0][1] = r2[3]
    
    v1[1] = v1[0] + h*f(s1[0], v1[0], 0)[2:]/2
    v2[1] = v2[0] + h*f(s2[0], v2[0], 0)[2:]/2
    
    M = int(len(tpoints-4)/2)
    
    for i in range(M):
        
        t = i*2
        
        r = s1[t+2] - s2[t+2]
        s1[t+2] = s1[t] + h*v1[t+1]
        k = h*f(r, v1[t+2], tpoints[t]+h)[2:]
        v1[t+2] = v1[t+1] + k/2
        v1[t+3] = v1[t+1] + k
        
        r = s2[t+2] - s1[t+2]
        s2[t+2] = s2[t] + h*v2[t+1]
        k = h*f(r, v2[t+2], tpoints[t]+h)[2:]
        v2[t+2] = v2[t+1] + k/2
        v2[t+3] = v2[t+1] + k
        
    return s1, s2

#initial conditions of each part r*=[[i],[ii],[iii]] ,[i]=[x,y] 
r1_all = np.array([[4.0, 4.0],[4.5, 4.0],[2.0, 3.0]])
r2_all = np.array([[5.2, 4.0],[5.2, 4.0],[3.5, 4.4]])

#ode solve parameters
a = 0.0
b = 1.0
N = 100
h = (b-a)/N

#solve and plot for each part
for i in range(len(r1_all)):
    #initial conditions & containers
    x1 = r1_all[i,0]
    y1 = r1_all[i,1]
    x2 = r2_all[i,0]
    y2 = r2_all[i,1]
    vx1 = 0.0
    vy1 = 0.0
    vx2 = 0.0
    vy2 = 0.0
    
    r1 = np.array([x1, y1, vx1, vy1])
    r2 = np.array([x2, y2, vx2, vy2])

    tpoints = np.arange(a,b,0.5*h)
    
    #solve ode
    s1, s2 = solve(r1, r2, func, h, tpoints)

    #extract x&y values from solution container
    x1 = s1[:,0]
    y1 = s1[:,1]
    x2 = s2[:,0]
    y2 = s2[:,1]
    
    #plot resulting function
    plt.figure()
    plt.plot(x1, y1, ".", label="1")
    plt.plot(x2, y2, ".", label="2")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Trajectories of molecules")
    plt.grid()
    plt.legend()
    path  = "../images/q2_b_" + "i"*i +".png"
    plt.savefig(path, dpi=600)




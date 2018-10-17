# -*- coding: utf-8 -*-
"""
Created on Tue Oct 16 18:37:32 2018

@author: Pierino
"""
#import libraries
import numpy as np
import pylab as plt

def func(A, B, t):
    """
    ODE function that computes the force on A due to B
    """
    #unpack coordinates of A & B into components 
    x1, y1, vx1, vy1 = A
    x2, y2, vx2, vy2 = B
    
    #compute relative coordinates
    x = (x2 - x1)
    y = (y2 - y1)
    
    #function constants
    epsilon = 1.0
    sigma = 1.0
    mass = 1.0

    #compute new values
    Dx = vx1
    Dy = vy2
    
    R = (x**2 + y**2)
    c1 = (48.0 * epsilon * R**(-7) * sigma**(12))/mass
    c2 = (-24.0* epsilon * R**(-4) * sigma**(6))/mass
    
    Dvx = (c1 + c2)*x + x1
    Dvy = (c1 + c2)*y + y1
    
    return np.array([Dx,Dy,Dvx,Dvy])

def plot(r1, r2, path, save):
    """
    Function that takes two 2D arrays with format [x,y,vx,vy] and plots the
    data sets on same graph.
    """
    #unpack data sets
    x1 = r1[:,0]
    y1 = r1[:,1]
    x2 = r2[:,0]
    y2 = r2[:,1]
    
    #plot both molecules
    plt.figure()
    plt.plot(x1, y1, ".", label="1")
    plt.plot(x2, y2, ".", label="2")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Trajectories of molecules")
    plt.grid()
    plt.legend()
    if save:
        plt.savefig(path, dpi=600)
        
def verlet(r1_initial ,r2_initial, f, h, tpoints):
    """
    Function that takes two 1D arrays, of the form [x,y,vx,vy], corresponding
    to the initial conditions of particles 1 and 2. Computes the integral
    of f using verlet method and returns two 2D arrays containting a set 
    [x, y, vx, vy] for each point.  
    """
    #set containers
    N = len(tpoints)
    r1 = np.zeros([N, 4])
    r2 = np.zeros([N, 4])
    
    #set initial conditions
    r1[0] = r1_initial
    r2[0] = r2_initial
    
    #compute first set
    r1[1,2:] = r1[0,2:] + h*f(r1[0], r2[0], t=0)[2:]/2
    r2[1,2:] = r2[0,2:] + h*f(r2[0], r1[0], t=0)[2:]/2
    
    
    # t=i => t + (h*n)/2=i+n  [t_step=0.005; h=0.01]
    
    for i in range(len(tpoints)//2-1):
        i = i*2
        r1[i+2,:2] = r1[i,:2] + h*r1[i+1,2:]
        k = h*f(r1[i+2], r2[i+2], tpoints[i+2])[2:]
        r1[i+2,2:] = r1[i+1,2:] + k/2
        r1[i+3,2:] = r1[i+1,2:] + k
        
        r2[i+2,:2] = r2[i,:2] + h*r2[i+1,2:]
        k = h*f(r2[i+2], r1[i+2], tpoints[i+2])[2:]
        r2[i+2,2:] = r2[i+1,2:] + k/2
        r2[i+3,2:] = r2[i+1,2:] + k
    
    s1 = []
    s2 = []    
    
    for i in range(len(r2)):
        if i%2 == 0:
           s1.append(r1[i])
           s2.append(r2[i])
    
    return np.array(s1), np.array(s2)    
        

#integration parameters
a = 0.0
b = 1.0
N = 100
h = (b-a)/N
tpoints = np.arange(a, b, 0.5*h)

#save images
save = False

#I
r1 = [4.0, 4.0, 0.0, 0.0]
r2 = [5.2, 4.0, 0.0, 0.0]
s1, s2 = verlet(r1, r2, func, h, tpoints)
path = "../images/q2_i.png"
plot(s1, s2, path, save)

#II
r1 = [4.5, 4.0, 0.0, 0.0]
r2 = [5.2, 4.0, 0.0, 0.0]
s1, s2 = verlet(r1, r2, func, h, tpoints)
path = "../images/q2_ii.png"
plot(s1, s2, path, save)

#III
r1 = [2.0, 3.0, 0.0, 0.0]
r2 = [3.5, 4.4, 0.0, 0.0]
s1, s2 = verlet(r1, r2, func, h, tpoints)
path = "../images/q2_iii.png"
plot(s1, s2, path, save)



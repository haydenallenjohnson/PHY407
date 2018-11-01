#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 31 10:09:32 2018

@author: Hayden

Program which uses the Gauss-Seidel method with overrelaxation to solve the
heat equation in steady-state on a rectangular domain (with a cut-out) given fixed boundary 
conditions. The program "animates" the solving process by plotting the values
at each iteration, before eventually plotting the final solution.
"""

import numpy as np
import matplotlib.pyplot as plt

#specify grid size and resolution
height = 8
width = 20
dx = 0.1

#specify target accuracy
target = 1e-6
w = 0.9

#create array to hold values
temp = np.zeros((int(height/dx)+1,int(width/dx)+1),float)

#specify boundary conditions of array
temp[-1,:51] = np.linspace(0,5,51)
temp[50:,50] = np.linspace(7,5,31)
temp[50,50:151] = 7
temp[50:,150] = np.linspace(7,5,31)
temp[-1,150:] = np.linspace(5,0,51)
temp[:,-1] = np.linspace(10,0,81)
temp[0,:] = 10
temp[:,0] = np.linspace(10,0,81)

#initialize value of max change from one iteration to the next
delta = 1.0

#run loop while precision remains less than desired
while delta > target:
    
    #create copy of old values of array to compare against later
    temp_old = np.copy(temp)

    #iterate over the rows and columns of the temperature array
    for i in range(temp.shape[0]):
        for j in range(temp.shape[1]):
            
            #if the point is on the boundary, don't do anything
            if i==0 or i==temp.shape[0]-1 or j==0 or j==temp.shape[1]-1:
                continue
            elif 50<=j<=150 and i>=50:
                continue
            #otherwise, we are not on the boundary, so compute the new value
            else:
                temp[i,j] = ((1+w)/4)*(temp[i+1,j] + temp[i-1,j] + temp[i,j+1] + temp[i,j-1]) - w*temp[i,j]

    #calculate the maximum difference of the new values from the old ones
    delta = max(abs(temp-temp_old))
    print(delta)
    
    #plot the new values (to make a frame in our bootleg "animation")
    plt.clf()
    plt.imshow(temp)
    plt.pause(0.01)
    plt.show()

#plot the final values of the temperature
plt.imshow(temp)
plt.title('Solution of Heat Equation with $\omega$='+str(w))
plt.xlabel('x position (mm)')
plt.ylabel('y position (mm)')
plt.colorbar(label='Temperature ($^\circ$C)')
plt.savefig('../images/q1_c.png') #b_'+str(w)[0]+'p'+str(w)[2]+'.png')

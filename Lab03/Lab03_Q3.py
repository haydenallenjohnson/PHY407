# -*- coding: utf-8 -*-
"""
Created on Wed Sep 26 11:08:07 2018

@author: Hayden and Pierino

Program reads altitude data out of file, computes gradient of the
altitude, computes the intensity of incident light, and plots the
altitude and light intensity
"""

#import modules
import numpy as np
import matplotlib.pyplot as plt
import struct

#define the size of the array of values, N
N = 1201

#create empty array to store height values
w = np.zeros((N,N))

#open file
f = open('N19W156.hgt', 'rb')

#iterate over the rows of the array
for i in range(N):
    
    #create empty array to store values for the row
    row = np.zeros(N)
    
    #iterate over the elements of the row
    for j in range(N):
        
        #read each value out of the file and store it to the row array
        buf = f.read(2)
        row[j] = struct.unpack('>h', buf)[0]
    
    #write the completed row into the 2d array storing the values
    w[i] = row

#close the file
f.close

#define function to compute gradient
def grad(w):
    
    #define value of h which specifies distance between grid points
    h = 90 #30m per arcsecond, 3 arcseconds per tile
    
    #create empty array to store gradient values
    v = np.zeros((N,N,2))
    
    #iterate over rows of the coordinate grid
    for i in range(N):
        
        #define an empty array to store gradient values for the row
        row = np.zeros((N,2))
        
        #iterate over the columns of the grid
        for j in range(N):

            #if in first column, compute dw/dx with forward derivative
            if j == 0:
                dw_dx = (w[i][j+1]-w[i][j])/h
                
            #if in last column, compute dw/dx with backward derivative
            elif j == (N-1):
                dw_dx = (w[i][j]-w[i][j-1])/h
                
            #otherwise, compute dw/dx using central difference
            else:
                dw_dx = (w[i][j+1]-w[i][j-1])/(2*h)
                
            #if in first row, compute dw/dy with forward derivative
            if i == 0:
                dw_dy = (w[i+1][j]-w[i][j])/h
                
            #if in last row, compute dw/dy with backward derivative
            elif i == (N-1):
                dw_dy = (w[i][j]-w[i-1][j])/h
                
            #otherwise, compute dw/dy using central difference
            else:
                dw_dy = (w[i+1][j]-w[i-1][j])/(2*h)
            
            #store computed partial derivatives to row array
            d_w = np.array([dw_dx,dw_dy])
            row[j]=d_w
            
        #store computed row to v array
        v[i] = row
        
    return v

#compute gradient
v = grad(w)

#define function to compute intensity
def I(grad_w, a_x, a_y, a_z):
    
    #extract values of dw/dx and dy/dx from gradient array
    dw_dx = grad_w[:,:,0]
    dw_dy = grad_w[:,:,1]
    
    #compute numerator and denominator
    numerator = a_x*dw_dx + a_y*dw_dy - a_z
    denominator = np.sqrt(dw_dx*dw_dx + dw_dy*dw_dy + 1)
    
    #compute value
    val = numerator/denominator
    
    #clip any negative results to zero
    val[val<0] = 0
    
    #return nalues
    return val

#compute light intensity over coordinate grid
intensity = I(v, -np.cos(np.pi), -np.sin(np.pi), 0)
intensity_2 = I(v, -np.cos(2*np.pi), -np.sin(2*np.pi), 0)

#specify max and min values of x and y in the grid
xmin = 155
xmax = 156
ymin = 19
ymax = 20

#create arrays of grid coordinates in x and y
x = np.linspace(xmax,xmin,1201)
y = np.linspace(ymax,ymin,1201)

#plot color map of altitude array using pcolormesh
plt.figure(0)
plt.pcolormesh(x, y, w, vmin=0)
plt.xlim([xmax,xmin])
plt.ylim([ymin,ymax])
plt.gca().set_aspect('equal', adjustable='box')
plt.suptitle('Altitude above sea level on Hawaii')
plt.xlabel('Longitude ($^\circ$W)')
plt.ylabel('Latitude ($^\circ$N)')
colorbar = plt.colorbar()
colorbar.set_label('Elevation (m)')
plt.savefig('altitude.png')

#plot color map of illumination intensity using pcolormesh
plt.figure(1)
plt.pcolormesh(x, y, intensity, vmax=0.2, cmap='gray')
plt.xlim([xmax,xmin])
plt.ylim([ymin,ymax])
plt.gca().set_aspect('equal', adjustable='box')
plt.title('Illumination of Hawaii at sunset')
plt.xlabel('Longitude ($^\circ$W)')
plt.ylabel('Latitude ($^\circ$N)')
plt.savefig('illumination.png')

#plot color map of illumination intensity using pcolormesh
plt.figure(2)
plt.pcolormesh(x, y, intensity_2, vmax=0.2, cmap='gray')
plt.xlim([xmax,xmin])
plt.ylim([ymin,ymax])
plt.gca().set_aspect('equal', adjustable='box')
plt.title('Illumination of Hawaii at sunrise')
plt.xlabel('Longitude ($^\circ$W)')
plt.ylabel('Latitude ($^\circ$N)')
plt.savefig('illumination_2.png')
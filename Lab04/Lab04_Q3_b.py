#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  3 11:43:38 2018

@author: Hayden

Program which, using both the regular relaxation and overrelaxation
methods, iterates over values of x starting from x=1 by computing f(x)
and taking this as the new x (or possibly with some adjustments, in
the case of overrelaxation). The program also prints the number of
iterations required to reach the specified accuracy for each of the
methods used.
"""

#import modules
import numpy as np

#define value of c, w for overrelaxation, and desired accuracy
c = 2.0
w = 0.7
acc = 1e-6

#define the function f
def f(x, c):
    return 1 - np.exp(-c*x)

#define a function for the deriviative of f
def f_1(x, c):
    return c*np.exp(-c*x)

#define a function to compute the error for the overrelaxation method
#note that taking w=0 corresponds to regular relaxation
def error(x1,x2,c,w):
    num = x1-x2
    a = (1+w)*f_1(x1,c) - w
    
    #if the value a is 0, the denominator will be very large, making
    #the error very small, so return an error small enough to halt the
    #iteration
    if a==0:
        return 1e-8
    
    #if the value of a is 1, the denominator will be very small, making 
    #the error very large, so return a large error that will not halt
    #the iteration
    elif a==1:
        return 1

    #otherwise, compute the value of the denominator as usual
    else:
        denom = 1 - (1/a)
        return abs(num/denom)

#------------------------------------------------------------------
#iterate over x using regular relaxation

#set initial value to iterate from, and compute new value
x_old = 1
x_new = f(x_old, c)

#define an iteration counter
i = 1

#while the error remains greater than desired accuracy, replace the
#old old value by the old new value and calculate the new new value 
#from the old new value; incrase iteration count by 1
while error(x_old, x_new, c, 0) > acc:
    x_old = x_new
    x_new = f(x_new, c)
    i += 1

#print iteration count
print('Relaxation solution: '+str(f(x_new,c)))
print('Number of iterations: '+str(i))


#------------------------------------------------------------------
#iterate over x using overrelaxation

#set initial value to iterate from, and compute new value
x_old = 1
x_new = f(x_old, c)

#define an iteration counter
j = 1

#while the error remains greater than desired accuracy, replace the
#old old value by the old new value and calculate the new new value 
#from the old new value; increase iteration count by 1
while error(x_old, x_new, c, w) > acc:
    x_old = x_new
    x_new = (1 + w)*f(x_new, c) - w*x_new
    j += 1

#print iteration count
print('Relaxation solution: '+str(f(x_new,c)))
print('Number of iterations: '+str(j))

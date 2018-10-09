#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 14:11:32 2018

@author: Hayden

Program which computes an estimate of the non-trivial root of 5e^{-x} - x + 5
using the binary search method, the (over)relaxation method, and Newton's
method to a desired accuracy, and also records the number of iterations
required to reach the result to the desired accuracy.
"""

#import modules
import numpy as np

#define the function f which we are solving for binary search method
#and for newton's method
def f(x):
    return 5*np.exp(-x) + x - 5

#define derivative of f to be used in newton's method
def f_1(x):
    return -5*np.exp(-x) + 1

#define the function g to use for solving with relaxation method
def g(x):
    return 5 - 5*np.exp(-x)

#define derivative of function for relaxation method
def g_1(x):
    return 5*np.exp(-x)

#define error for the (over)relaxation method
def relax_error(x1,x2,w):
    num = x1-x2
    a = (1+w)*g_1(x1) - w
    
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

#define error for newton's method
def newton_error(x_old,x_new):
    return abs(x_new - x_old)

#define desired accuracy
acc = 1e-6

#-------------------------------------------------------------
#sove using binary search method

#set initial interval
x1 = 4
x2 = 6

#confirm that we have not chosen an interval which does contain a root
if np.sign(f(x1)) == np.sign(f(x2)):
    print("Warning: invalid interval")

#compute f at the boundaries
f_x1 = f(x1)
f_x2 = f(x2)

#initiate iteration counter
i = 0

#while the distance between x1 and x2 is greater than the desired accuracy:
while x2-x1 > 1e-6:
    
    #increase iteration count
    i += 1
    
    #compute the x value of the midpoint of the current interval
    x_mid = (x1+x2)/2
    
    #if f(x_mid) has same sign as f(x1), make x_mid the new x1
    if np.sign(f(x_mid)) == np.sign(f(x1)):
        x1 = x_mid
    
    #if f(x_mid) has same sign as f(x2, make x_mid the new x2)
    elif np.sign(f(x_mid)) == np.sign(f(x2)):
        x2 = x_mid
    
    #if we happend to get lucky and hit the solution, terminate the loop
    elif np.sign(f(x_mid)) == 0:
        break

#calculate the midpoint one last time and take this as the final estimate
#of the solution
x_mid = (x1+x2)/2
i+=1

#print results
print('Binary solution: '+str(x_mid))
print('Number of iterations: '+str(i))

#-------------------------------------------------------------
#sove using relaxation method

#set value of w
w = 0.04

#set initial value to iterate from, and compute new value
x_old = 4
x_new = g(x_old)

#define an iteration counter
j = 1

#while the error remains greater than desired accuracy, replace the
#old old value by the old new value and calculate the new new value 
#from the old new value; incrase iteration count by 1
while relax_error(x_old, x_new, w) > acc:
    x_old = x_new
    x_new = (1 + w)*g(x_new) - w*x_new
    j += 1

#print results
print('Relaxation solution: '+str(g(x_new)))
print('Number of iterations: '+str(j))

#-------------------------------------------------------------
#sove using Newton's method

# define initial guess for x, and compute new value
x_old = 6
x_new = x_old - (f(x_old)/f_1(x_old))

#define iteration counter
k = 1

#while error on x_old remains greater than desired accuracy, compute
#new values of x
while newton_error(x_old,x_new) > acc:
    k += 1
    x_old = x_new
    x_new = x_old - (f(x_old)/f_1(x_old))
    
#print results
print('Newton solution: '+str(x_new))
print('Number of iterations: '+str(k))

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 12 21:44:52 2018

@author: Hayden
"""

#import modules
import numpy as np
import matplotlib.pyplot as plt

#define function to multiply two matrices together
def multiply(A, B):
    
    #check to make sure that the arrays are the right sizes to be multiplied
    if len(A[0]) != len(B):
        return float('nan')
    
    #if the test is passed, carry on
    else:
        #define empty result matrix
        C = np.empty((len(A),len(B[0])))
        
        #iterate over the rows of A
        for i in range(len(A)):
            #iterate over the columns of B
            for j in range(len(B[0])):
                #(re)initialize the variable to store the sum
                s=0
                #iterate over the length of the row of A (column of B)
                for k in range(len(A[0])):
                    #add the product of the terms to the sum
                    s += A[i,k]*B[k,j]
                    
                #store the value to the array
                C[i,j] = s
                
        #return the array
        return C

# -*- coding: utf-8 -*-
"""
Computes the standard deviation of a given data set using a one pass and two
pass method and compares their relative error values, using a numpy function
to compute the 'true' value.

Created on Tue Sep 18 14:07:16 2018
@author: Pierino Zindel
"""

#Import libraries
import numpy as np


#Define standard deviation function that computes mean first (i.e. two pass)
def stdev_two_pass_method(data_points):
    #Compute the mean
    N = len(data_points)
    data_total = 0
    for x in data_points:
        data_total = data_total + x
    data_ave = data_total / N
    
    #Compute and return the standard deviation
    sum_1 = 0
    for x in data_points:
        sum_1 = sum_1 + (x - data_ave)**2
    sigma = (sum_1 / (N-1))**0.5
    return sigma
    

#Define standard deviation function that uses one pass through data set
def stdev_one_pass_method(data_points):
    #Compute the sum
    N = len(data_points)
    sum_1 = 0
    data_total = 0
    for x in data_points:
        data_total = data_total + x
        sum_1 = sum_1 + x**2
    
    #Compute the standard deviation
    data_ave = data_total/N
    var = (sum_1 - N*data_ave**2) / (N-1)
    if var < 0:
        raise Exception("Variance is negative")
    sigma = var**0.5
    return sigma
    
#Define function that returns the relative error between a given standard
#deviation and the 'true' standard deviation
def stdev_rel_error(data_points, sigma_est):
    sigma_true = np.std(data_points, ddof=1)
    return abs((sigma_est - sigma_true) / sigma_true)
    
#Open data set file
source_file = "cdata.txt"
speed_c = np.loadtxt(source_file, usecols=(0), 
                         unpack=True)

#Compute standard deviation using each method
sigma_two_pass = stdev_two_pass_method(speed_c)
sigma_one_pass = stdev_one_pass_method(speed_c)

#Compute relative error of each method
rel_error_two_pass = stdev_rel_error(speed_c, sigma_two_pass)
rel_error_one_pass = stdev_rel_error(speed_c, sigma_one_pass)

#Print results
print("The relative error of method/eqtn one is: ", rel_error_two_pass)
print("The relative error of method/eqtn two is: ", rel_error_one_pass)

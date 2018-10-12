#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 14:21:08 2018

@author: Hayden

Program which imports DOW closing value data from the dow.txt file
and plots it, and then takes the fourier transform of the data, sets
the last 90% (and then 98%) of the fourier coefficients to be zero,
and takes the inverse fourier transform to produce smoothed functions
representing only the low-frequency components of the data.
"""

#import modules
import numpy as np
import matplotlib.pyplot as plt

#------------------------------------------------------------------
#part a

#import closing value data from dow.txt
dow = np.loadtxt('../data/dow.txt')

#define an array of the days since the start of the data
days = np.arange(0,len(dow),1)

#plot dow closing value as a function of days
plt.figure(0)
plt.plot(days, dow)
plt.title('DOW closing value over time')
plt.xlabel('Business days since start of the data stream')
plt.ylabel('DOW closing value (in points)')
plt.savefig('../images/1b_dow.png')

#------------------------------------------------------------------
#part b

#compute fourier coefficients of discrete fourier transform
c = np.fft.rfft(dow)

#------------------------------------------------------------------
#part c

#set the fraction of the fft coefficients that we want to keep non-zero
keep_frac = 0.1

#calculate the number of values of c to keep non-zero
keep_length = round(keep_frac*len(c))

#create a new array from c with only the first designated fraction 
#non-zero
c_new = np.zeros(len(c), dtype=complex)
c_new[:keep_length] = c[:keep_length]

#------------------------------------------------------------------
#part d

#compute inverse fourier transform of filtered c values
dow_filtered = np.fft.irfft(c_new)

#plot both the original and filtered data
plt.figure(1)
plt.plot(days, dow, label='Original')
plt.plot(days, dow_filtered, label='Filtered (10%)')
plt.legend(loc='lower left')
plt.title('DOW closing value with 10% filter')
plt.xlabel('Business days since start of the data stream')
plt.ylabel('DOW closing value (in points)')
plt.savefig('../images/1b_filtered_10.png')

#------------------------------------------------------------------
#part e

#set the fraction of the fft coefficients that we want to keep non-zero
keep_frac_2 = 0.02

#calculate the number of values of c to keep non-zero
keep_length_2 = round(keep_frac_2*len(c))

#create a new array from c with only the first designated fraction 
#non-zero
c_new_2 = np.zeros(len(c), dtype=complex)
c_new_2[:keep_length_2] = c[:keep_length_2]

#compute inverse fourier transform of filtered c values
dow_filtered_2 = np.fft.irfft(c_new_2)

#plot both the original and filtered data
plt.figure(2)
plt.plot(days, dow, label='Original')
plt.plot(days, dow_filtered_2, label='Filtered (2%)')
plt.legend(loc='lower left')
plt.title('DOW closing value with 2% filter')
plt.xlabel('Business days since start of the data stream')
plt.ylabel('DOW closing value (in points)')
plt.savefig('../images/1b_filtered_02.png')

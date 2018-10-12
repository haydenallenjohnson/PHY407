#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 21:18:10 2018

@author: Hayden
"""

#import modules
import numpy as np
import matplotlib.pyplot as plt
from dcst import dct, idct

#import closing value data from dow.txt
dow = np.loadtxt('../data/dow2.txt')

#define an array of the days since the start of the data
days = np.arange(0,len(dow),1)

#compute fourier coefficients of discrete fourier transform
c = np.fft.rfft(dow)

#set the fraction of the fft coefficients that we want to keep non-zero
keep_frac = 0.02

#calculate the number of values of c to keep non-zero
keep_length = round(keep_frac*len(c))

#create a new array from c with only the first designated fraction 
#non-zero
c_new = np.zeros(len(c), dtype=complex)
c_new[:keep_length] = c[:keep_length]

#compute inverse fourier transform of filtered c values
dow_filtered = np.fft.irfft(c_new)

#plot both the original and filtered data
plt.figure(0)
plt.plot(days, dow, label='Original')
plt.plot(days, dow_filtered, label='Filtered (DFT)')
plt.legend(loc='upper left')
plt.title('DOW closing value raw and filtered')
plt.xlabel('Business days since start of the data stream')
plt.ylabel('DOW closing value (in points)')
plt.savefig('../images/1c_dft.png')

#-------------------------------------------------------------------
#part b

#compute dct
a = dct(dow)

#calculate the number of values of a to keep non-zero
keep_length_2 = round(keep_frac*len(a))

#create a new array from c with only the first designated fraction 
#non-zero
a_new = np.zeros(len(a), dtype=complex)
a_new[:keep_length_2] = a[:keep_length_2]

#compute inverse fourier transform of filtered c values
dow_filtered_2 = idct(a_new)

#plot both the original and filtered data
plt.figure(2)
plt.plot(days, dow, label='Original')
plt.plot(days, dow_filtered_2, label='Filtered (DCT)')
plt.legend(loc='upper left')
plt.title('DOW closing value with 2% filter')
plt.xlabel('Business days since start of the data stream')
plt.ylabel('DOW closing value (in points)')
plt.savefig('../images/1c_dct.png')

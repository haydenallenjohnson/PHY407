# -*- coding: utf-8 -*-
"""
Takes a data set corresponding to a blurred image and deconvolves it using 
fourier transforms and a Gaussian point spread function to produce a clearer
image.

Created on Wed Oct 10 09:38:31 2018
@author: Pierino Zindel
"""

#import libraries and functions
import numpy as np
import matplotlib.pyplot as plt
from numpy.fft import rfft2, irfft2

#define Gaussian point spread function
def Gaussian(x,y,sigma):
    return np.exp(-(x**2 + y**2)/(2.0*sigma**2))

#read in blurred data file
source_file = "blur.txt"
with open(source_file) as file:
    data = np.array([[float(digit) for digit in line.split()] for line in file])
    
#create density plot of blurred data
plt.figure()
plt.imshow(data, cmap='gray')
bar = plt.colorbar()
bar.set_label("greyscale")
plt.title("Blurry image")
plt.xlabel("x")
plt.ylabel("y")
plt.savefig("../images/blurred_image.png", dpi=600)
plt.show()

#define variables for spread function computation
N = len(data)
M = len(data[0])
gauss_grid = np.empty([N,M])
sigma = 25

#compute spread function for each point in the image
for i in range(N):
    x = i
    if x > N/2:
        x -= N
    for j in range(M):
        y = j
        if y > M/2:
            y -= M        
        gauss_grid[i,j] = Gaussian(x,y,sigma)

#create density plot of spread function
plt.figure()
plt.imshow(gauss_grid, cmap='gray')
bar = plt.colorbar()
bar.set_label("brightness")
plt.title("Gaussian point spread function")
plt.xlabel("x")
plt.ylabel("y")
plt.savefig("../images/spread_function.png", dpi=600)
plt.show()

#fourier transform the data and spread function
data_fft = rfft2(data)
gauss_fft = rfft2(gauss_grid)

#define variables for computed fft of deconvolved image
J = len(gauss_fft)
K = len(gauss_fft[0])
epsilon = 1e-3
sharp_data_fft = np.zeros((J,K), dtype=np.complex_)

#apply spread function to data to retrieve fft of deconvolved image
for i in range(J):
    for j in range(K):
        if gauss_fft[i,j] > epsilon:
            sharp_data_fft[i,j] = data_fft[i,j]/gauss_fft[i,j]
        else:
            sharp_data_fft[i,j] = data_fft[i,j]

#apply inverse fourier transform to retrieve clear image data
sharp_data = irfft2(sharp_data_fft)

#create density plot of deconvolved image
plt.figure()
plt.imshow(sharp_data, cmap='gray')
bar = plt.colorbar()
bar.set_label("greyscale")
plt.title("Deconvolved imgae")
plt.xlabel("x")
plt.ylabel("y")
plt.savefig("../images/deconvolved_image.png", dpi=600)
plt.show()

# -*- coding: utf-8 -*-
"""
Computes the wavefunctions of the quantum harmonic oscillator using Hermite 
polynomials and then produces the energy of the wavefunctions and the 
epxectation values and uncertainties of the position-squared
and momentum-squared via Gaussian quadrature integration. 

Created on Wed Sep 26 09:26:56 2018
@author: Pierino Zindel
"""

#import libraries
import numpy as np
import pylab as plt
from gaussxw import gaussxwab
"""
NOTE: The gaussxw module was written by Mark Newmann for  
'Computational Physics', and downloaded from the site: 
http://www-personal.umich.edu/~mejn/cp/programs/gaussxw.py 
"""

#computes hermite polynomial of order n for a given position
def H(n,x):
    if n == -1:
        return 0
    elif n == 0:
        return 1
    elif n < 0:
        raise Exception("Hermite polynomial undefined for n<0")
    elif n == 1:
        return 2*x
    else:
        return 2*x*H(n-1,x) - 2*(n-1)*H(n-2,x)

#computes the coefficient for the wavefunction for a given n and position
def wavefunction_coeff(n,x):
    a = np.exp(-x**2/2)
    b = (float(2**n)*float(np.math.factorial(n))*(np.pi**0.5))**0.5
    return a/b

#computes the wavefunction for a given n and position
def wavefunction(n,x):
    return wavefunction_coeff(n,x) * H(n,x)

#computes the derivative of the wavefunction for a given n and position
def wavefunc_derivative(n,x):
    return wavefunction_coeff(n,x) * (-x*H(n,x) + 2*n*H(n-1,x))

#computes the expectation value of the position square for a given n
def expect_x2(n, p):
    #function to be integrated over
    def func(x):
        return ((1+x**2)/((1-x**2)**2))*((x/(1-x**2))**2)*(np.abs(wavefunction(n,(x/(1-x**2)))))**2
    N,a,b = p
    I = 0.0
    x,w = gaussxwab(N,a,b)
    for i in range(N):
        I += w[i]*func(x[i])
    return I

#computes the expectation value of the momentum square for a given n
def expect_p2(n, p):
    #function to be integrated over
    def func(x):
        return ((1+x**2)/((1-x**2)**2))*(np.abs(wavefunc_derivative(n,(x/(1-x**2)))))**2
    N,a,b = p
    x,w = gaussxwab(N,a,b)
    I = 0.0
    for i in range(N):
        I += w[i]*func(x[i])
    return I

#computes the energy given the expectation values of the position and momentum square
def energy(x2, p2):
    return 0.5*(x2+p2)

#computes the uncertainty of position for a given n
def x_uncer(n, p):
    return np.sqrt(expect_x2(n, p))
    
#computes the uncertainty of momentum for a given n
def p_uncer(n, p):
    return np.sqrt(expect_p2(n, p))

#compute the wavefunctions for n=0,1,2,3 and plot them
x = np.linspace(-4, 4, 1000)
plt.figure()
plt.grid()
plt.ylabel("$\Psi_n$(x)")
plt.xlabel("x")
plt.title("Harmonic Oscillator Wavefunctions")
for i in range(4):
    y = wavefunction(i,x)
    wave_label = "n=" + str(i)
    plt.plot(x, y, label=wave_label)
plt.legend()
plt.savefig("wavefunc_n4.png")

#compute the wavefunction for n=30 and plot it
x = np.linspace(-10, 10, 1000)
y = wavefunction(30,x)
plt.figure()
plt.grid()
plt.ylabel("$\Psi_{30}$(x)")
plt.xlabel("x")
plt.title("Harmonic Oscillator Wavefunction for n=30")
plt.plot(x,y,label="n=30")
plt.legend()
plt.savefig("wavefunc_n30.png")

#parameters for integration using gaussian method [N,a,b], where N is the
#number of points, a and b are the limits [Note: limits of +/-1 are 
#equivalent to +/- infinity]
p = [100, -1.0, 1.0]

#compute the energy and uncertainty in position and momentum for n=0,...,15
for i in range(16):
    x2 = expect_x2(i,p)
    p2 = expect_p2(i,p)
    print("\nUncertainty for n="+str(i)+" is: ")
    print("\nx: ", np.sqrt(x2))
    print("\np: ", np.sqrt(p2))
    print("\nEnergy for n="+str(i)+" is: ", energy(x2,p2),"\n")
    


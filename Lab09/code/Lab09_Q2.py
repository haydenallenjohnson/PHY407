#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 14 10:16:44 2018

@author: Hayden
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy.constants as pc

#define path to ffmpeg for animation
plt.rcParams['animation.ffmpeg_path'] = '/Users/Hayden/ffmpeg-4.1-macos64-static/bin/ffmpeg'
#zeranoe

#define constants
m = pc.m_e #electron mass
hbar = pc.hbar
e = pc.e #elementary charge
eps0 = pc.epsilon_0
pi = np.pi

#define number of samples and stepsize
N = 1000
L = 1e-8
a = L/N

#define timestep and final time
timesteps = 5000
h = 1e-18
t_final = timesteps*h

#define array of interior x values
x = np.arange(0, L+a, a)

#define values
x0 = L/2
sigma = 1e-10
kappa = 5e10

#calculate matrix element values
a1 = np.complex(1, h*hbar/(2*m*a*a))
a2 = np.complex(0, -h*hbar/(4*m*a*a))

b1 = np.complex(1, -h*hbar/(2*m*a*a))
b2 = np.complex(0, h*hbar/(4*m*a*a))

#constuct matrices A and B
A = np.zeros((N-1,N-1), dtype='complex')
B = np.zeros((N-1,N-1), dtype='complex')
for i in range(N-1):
    for j in range(N-1):
        if j==i:
            A[i,j] = a1
            B[i,j] = b1
        elif j==i-1 or j==i+1:
            A[i,j] = a2
            B[i,j] = b2

#compute inverse of A
Ainv = np.linalg.inv(A)

#compute inverse of A dotted with B
AinvB = np.dot(Ainv, B)

#compute psi at t=0
psi0 = np.exp(-((x-x0)**2)/(2*sigma*sigma))*np.exp(1j*kappa*x)

#specify boundary conditions
psi0[0] = 0
psi0[-1] = 0

#create empty array to store wavefunction values as function of time
wavefunction = np.zeros((timesteps,len(x)), dtype='complex')
wavefunction[0] = psi0

#iterate over timesteps and compute new value of psi
for i in range(timesteps-1):
    psi_old = wavefunction[i,1:-1]
    psi_new = np.dot(AinvB, psi_old)
    wavefunction[i+1,1:-1] = psi_new

#create axes for animation
fig = plt.figure(0)

#specify plots for animation (at every fifth timestep)
plots = []
for i in range(timesteps):
    if i%5==0:
        plots.append(plt.plot(x, np.real(wavefunction[i]), 'b'))

ani = animation.ArtistAnimation(fig, plots, interval=2, blit=True)
#save the animation
FFwriter = animation.FFMpegWriter(fps=30)
ani.save('../images/animation.mp4', writer = FFwriter)

plt.figure(1)
plt.plot(x,np.real(psi0))
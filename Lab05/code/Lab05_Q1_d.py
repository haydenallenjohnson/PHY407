# -*- coding: utf-8 -*-
"""
Computes the Fourier transform of waveforms corresponding to a note played on 
a piano and a trumpet recorded at a sampling rate of 44.1kHz and then outputs 
the frequency of the note.

Created on Thu Oct 11 16:40:40 2018
@author: Pierino Zindel
"""

#import libraries
import numpy as np
import pylab as plt

#read in data files
source_file = "../data/piano.txt"
source_file_2 = "../data/trumpet.txt"
p_data = np.loadtxt(source_file, unpack=True)
t_data = np.loadtxt(source_file_2, unpack=True)

#define ranges for the data sets
N = len(p_data)
M = len(t_data)
x = np.linspace(1, N, N)
z = np.linspace(1, M, M)

#plot data for piano waveform
plt.figure()
plt.title("Paino Data")
plt.xlabel("x")
plt.ylabel("y")
plt.grid()
plt.plot(x, p_data, label='paino waveform')
plt.legend()
plt.savefig("../images/piano_wave.png", dpi=600)

#plot data for trumpet waveform
plt.figure()
plt.title("Trumpet Data")
plt.xlabel("x")
plt.ylabel("y")
plt.grid()
plt.plot(z, t_data, label='trumpet waveform')
plt.legend()
plt.savefig("../images/trumpet_wave.png", dpi=600)

#perfrom fft on both data sets
p_fft = np.fft.fft(p_data)
t_fft = np.fft.fft(t_data)

#define a new range for plotting
N = 10000
x = np.linspace(1, N, N)

#plot data for fft of piano waveform
plt.figure()
plt.title("Paino Fouier Transform for First $1e4$ Coefficients")
plt.xlabel("Fourier coefficient number (k)")
plt.ylabel("$|c_k|$")
plt.grid()
plt.plot(x, abs(p_fft[:N]), label='fft of piano waveform')
plt.legend()
plt.savefig("../images/piano_fft.png", dpi=600)

#plot data for fft of trumpet waveform
plt.figure()
plt.title("Trumpet Fouier Transform for First $1e4$ Coefficients")
plt.xlabel("Fourier coefficient number (k)")
plt.ylabel("$|c_k|$")
plt.grid()
plt.plot(x, abs(t_fft[:N]), label='fft of trumpet waveform')
plt.legend()
plt.savefig("../images/trumpet_fft.png", dpi=600)

#number of bins for our fft
p_bin_num = len(p_data)/2
t_bin_num = len(p_data)/2

#sampling rate of fft of recording
fft_sample_rate = 44100/2

#frequency resolution per fft bin
p_freq_res = fft_sample_rate/p_bin_num
t_freq_res = fft_sample_rate/t_bin_num

#multiplying our fft coefficient numbers by our fft freq resolution should 
#yield the original frequency

#find frequencies of largest peaks
p_freq = p_freq_res*x[np.argmax(p_fft[:N])]
t_freq = p_freq_res*x[np.argmax(t_fft[:N])]

#print frequencies
print("\nThe frequency of the piano note is: \n", p_freq)
print("\nThe frequency of the trumpet note is: \n", t_freq)

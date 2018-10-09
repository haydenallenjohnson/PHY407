# -*- coding: utf-8 -*-
"""
Solves a linear system corresponding to a RC and RLC circuit to compute the 
voltages of the circuit at three different points. Utilizes a partial 
pivoting algorithm to solve the linear system.

Created on Thu Oct  4 17:06:17 2018
@author: Pierino Zindel
"""

#import libraries
import numpy as np
import SolveLinear as sl
import pylab as plt

#define constants
R1 = 1000 #ohms
R3 = 1000 #ohms
R5 = 1000 #ohms
R2 = 2000 #ohms
R4 = 2000 #ohms
R6 = 2000 #ohms
C1 = 1*10**(-6) #farad
C2 = 0.5*10**(-6) #farad
xplus = 3 #volts
omega = 1000 #rads per second
L = R6/omega #henry     

#define matrix entries
a11 = complex(1.0/R1 + 1.0/R4, omega*C1)
a12 = complex(0, -omega*C1)
a13 = complex(0, 0)
a21 = complex(0, -omega*C1)
a22 = complex(1.0/R2 + 1.0/R5, omega*C1 + omega*C2)
a23 = complex(0, -omega*C2)
a31 = complex(0, 0)
a32 = complex(0, -omega*C2)
a33 = complex(1.0/R3 + 1.0/R6, omega*C2)
a33_ALT = complex(1.0/R3, -1.0/(omega*L) + omega*C2)
v1 = complex(xplus/R1, 0)
v2 = complex(xplus/R2, 0)
v3 = complex(xplus/R3, 0)

#create matrices
A = np.array([[a11, a12, a13], [a21, a22, a23], [a31, a32, a33]])
A_ALT = np.array([[a11, a12, a13], [a21, a22, a23], [a31, a32, a33_ALT]])
v = np.array([v1, v2, v3])

#call partial pivoting to solve linear systems
x = sl.PartialPivot(A, v)
x_ALT = sl.PartialPivot(A_ALT, v)

#define function to compute voltages
def voltage(x, t):
    V = []
    if isinstance(t,int):
        return x*np.exp(complex(0, omega*t))
    for i in range(len(t)):
        V.append(x*np.exp(complex(0, omega*t[i])))
    return V

#compute voltages at t=0
t = 0
V1 = voltage(x[0], t)
V2 = voltage(x[1], t)
V3 = voltage(x[2], t)
V1_ALT = voltage(x_ALT[0], t)
V2_ALT = voltage(x_ALT[1], t)
V3_ALT = voltage(x_ALT[2], t)

#print the voltages and phases at t=0 for the RC circuit
print("|V1| = ", abs(V1))
print("phase x1 = ", np.angle(x[0]))
print("\n|V2| = ", abs(V2))
print("phase x2 = ", np.angle(x[1]))
print("\n|V3| = ", abs(V3))
print("phase x3 = ", np.angle(x[2]))

#compute voltages for two cycles for the RC circuit
t = np.linspace(0, 0.012, 100)
V1 = voltage(x[0], t)
V2 = voltage(x[1], t)
V3 = voltage(x[2], t)

#plot real component of the voltages for the RC circuit
plt.figure()
plt.title("RC Cicuit Voltages")
plt.xlabel("Time, t(s)")
plt.ylabel("Voltage, V(V)")
plt.grid()
plt.plot(t, np.real(V1), label="$V_1=x_1e^{i\omega t}$", color='blue')
plt.plot(t, np.real(V2), label="$V_2=x_2e^{i\omega t}$", color='red')
plt.plot(t, np.real(V3), label="$V_3=x_3e^{i\omega t}$", color='orange')
plt.legend()
plt.savefig("RC.png")

#print the voltages and phases at t=0 for RLC circuit
print("|V1_RLC| = ", abs(V1_ALT))
print("phase x1_RLC = ", np.angle(x_ALT[0]))
print("\n|V2|_RLC = ", abs(V2_ALT))
print("phase x2_RLC = ", np.angle(x_ALT[1]))
print("\n|V3|_RLC = ", abs(V3_ALT))
print("phase x3_RLC = ", np.angle(x_ALT[2]))

#compute voltages for two cycles for the RLC circuit
t = np.linspace(0, 0.012, 100)
V1_ALT = voltage(x_ALT[0], t)
V2_ALT = voltage(x_ALT[1], t)
V3_ALT = voltage(x_ALT[2], t)

#plot the real component of the voltages for the RLC circuit
plt.figure()
plt.title("RLC Cicuit Voltages")
plt.xlabel("Time, t(s)")
plt.ylabel("Voltage, V(V)")
plt.grid()
plt.plot(t, np.real(V1_ALT), label="$V_1=x_1e^{i\omega t}$", color='blue')
plt.plot(t, np.real(V2_ALT), label="$V_2=x_2e^{i\omega t}$", color='red')
plt.plot(t, np.real(V3_ALT), label="$V_3=x_3e^{i\omega t}$", color='orange')
plt.legend()
plt.savefig("RLC.png")

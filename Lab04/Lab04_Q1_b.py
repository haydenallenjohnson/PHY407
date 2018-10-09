# -*- coding: utf-8 -*-
"""
Solves linear systems composed of random square matrices of sizes 5-200
using Gaussian elimination, LU decomposition, and partial pivoting; and
compares the computation times and errors of each method.

Created on Wed Oct  3 09:41:11 2018
@author: Pierino Zindel
"""

#import libraries
import numpy as np
from time import time
import SolveLinear as sl
import pylab as plt

#define variables
N = 1000
A = []
v = []
x = []
error_1 = []
t_duration_1 = []
error_2 = []
t_duration_2 = []
error_3 = []
t_duration_3 = []

#create series of matrices with random entries
for i in range(5, N, 10):
    A.append(np.random.rand(i,i))
    v.append(np.random.rand(i))

M = len(v)

#gaussian elimination execution and timer
for i in range(M):
    start = time()
    x = sl.GaussElim(A[i], v[i])
    stop = time()
    t_duration_1.append(stop-start)
    v_sol = np.dot(A[i],x)
    error_1.append(np.mean(np.abs(v[i] - v_sol)))

#partial pivoting execution and timer
for i in range(M):
    start = time()
    x = sl.PartialPivot(A[i], v[i])
    stop = time()
    t_duration_2.append(stop-start)
    v_sol = np.dot(A[i],x)
    error_2.append(np.mean(np.abs(v[i] - v_sol)))

#LU decomp execution and timer
for i in range(M):
    start = time()
    x = np.linalg.solve(A[i], v[i])
    stop = time()
    t_duration_3.append(stop-start)
    v_sol = np.dot(A[i],x)
    error_3.append(np.mean(np.abs(v[i] - v_sol)))

#plot the error and time data of each method on a log-log plot
plt.figure()
plt.title("Performance of Linear System Solution Methods")
plt.ylabel("Ave error of solution, $log(mean(|v_{computed}-v_{actual}|)$")
plt.xlabel("Computation time for solution, $log(t)$")
plt.grid()
plt.loglog()
plt.plot(t_duration_1, error_1, '.', label='Gaussian Elimination', color='blue')
plt.plot(t_duration_2, error_2, '.', label='Partial Pivoting', color='orange')
plt.plot(t_duration_3, error_3, '.', label='LU Decomposition', color='red')
plt.legend()
plt.savefig("lin_sol.png")
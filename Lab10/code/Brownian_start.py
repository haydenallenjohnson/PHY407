
#################################################################
# This program simulates Brownian motion in the presence of walls
# Note that the physical behaviour would be to stick to walls,
# which is the purpose of Q1a.
# Author: Nico Grisouard, University of Toronto
# Date: 14 November 2018
#################################################################
import numpy as np
import matplotlib.pyplot as plt


def nextmove(x, y):
    """ randomly choose a direction
    1 = up, 2 = down, 3 = left, 4 = right"""
    direction =  # COMPLETE

    if direction == 1:  # move up
        y += 1
    elif direction == 2:  # move down
        y -= 1
    elif direction == 3:  # move right
        x += 1
    elif direction == 4:  # move left
        x -= 1
    else:
        print("error: direction isn't 1-4")

    return x, y


# %% main program starts here ------------------------------------------------|
# YOU NEED TO FINISH IT!

plt.ion()

Lp = 101  # size of domain
Nt = 5000  # number of time steps
# arrays to record the trajectory of the particle
# COMPLETE

centre_point = (Lp-1)//2  # middle point of domain
xp = centre_point
yp = centre_point

for i range(Nt):
    xpp, ypp = nextmove(xp, yp)
    # AND OFF YOU GO!

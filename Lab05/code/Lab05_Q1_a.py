#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 16:25:26 2018

@author: Hayden
"""

import numpy as np
import matplotlib.pyplot as plt

x = np.arange(0,100,1)

y = x**2

plt.figure(0)
plt.plot(x,y)
plt.savefig('../images/name.png')
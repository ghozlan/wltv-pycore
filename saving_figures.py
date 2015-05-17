# -*- coding: utf-8 -*-
"""
Created on Sat May 16 12:47:39 2015

@author: Hassan
"""


import matplotlib.pylab as plt
import numpy as np

def plot():
    plt.subplot(2,1,1)
    plt.plot((np.arange(10))**2)
    plt.subplot(2,1,2)
    plt.plot((np.arange(10))*2)
    
plot()
plt.savefig('test2.png')
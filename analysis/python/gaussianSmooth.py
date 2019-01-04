#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 28 12:32:33 2018

@author: ty
"""
import numpy as np
import matplotlib.pyplot as plt

win = 0.66667 #width of smoothing window in THz

#first, read data in from tr.dat file using readTR.dat

thz = data[:,0]
dom = np.real((thz[1]-thz[0])*2*np.pi*1e12)

### GAUSSIAN SMOOTHING ######
gwin = round(win*1e12*2*np.pi/dom) #number of array elements in window
if gwin % 2 == 0: #make sure its odd sized array
    gwin = gwin+1
if gwin == 1:
    gauss = np.array([1]) #if array is size 1, convolve with self
else:
    n = 2*np.arange(0,gwin)/(gwin-1)-(1) #centered at 0, sigma = 1
    n = (3*n) #set width of profile to 6*sigma i.e. end values ~= 0
    gauss = np.exp(-np.multiply(n,n)/2.0) #gaussian fx to convolve with
    gauss = gauss/np.sum(gauss) #normalized gaussian

#

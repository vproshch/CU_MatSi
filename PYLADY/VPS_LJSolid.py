#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
Calculate vibrational power spectrum of an argon-like lennard jones solid.
"""

### VARIABLES ###
import numpy as np
import matplotlib.pyplot as plt
import ty

velsFile = 'vels.dat'
posFile = 'data.argon'

steps = 500000 #run time of MD simulations
dt = 5e-15 #MD time step
dn = 100 #print frequency
mass = 39.95 #argon
a = 5.376 #lattice constant of argon

split = 2 #times to split data for averaging
prints = steps/dn #times data is printed
tn = prints/split #timesteps per chunk

pi = np.pi #tired of forgetting the 'np' part...

om = np.arange(0,tn)*2*np.pi/(tn*dt*dn) #angular frequency
thz = om/2/np.pi*1e-12 #convert to THz
dom = om[1]-om[0]
win = 0.01 #gaussian smoothing window
####

### GET COORDS ###
with open(posFile, 'r') as fid:
    for i in range(2): #skip comments
            fid.readline()
    num = int(fid.readline().strip().split()[0]) #number of atoms
    fid.readline() #skip comments
    types = int(fid.readline().strip().split()[0])
#    for i in range(12): #skip comments
#        fid.readline()
#    pos = np.zeros((num,5))
#    for i in range(num):
#        tmp = fid.readline().strip().split()
#        pos[i,0] = int(tmp[0]) #atom id
#        pos[i,1] = int(tmp[1]) #atom type
#        pos[i,2] = float(tmp[2]) #x-coordinate
#        pos[i,3] = float(tmp[3]) #y-coordinate
#        pos[i,4] = float(tmp[4]) #z-coordinate
#del posFile
####

## GET VELOCITIES ###
velsFFT = np.zeros((tn,num,3)).astype(complex) #FFT of velocity data
with open(velsFile, 'r') as fid: 
    vels = np.zeros((tn,num,3))
    for i in range(split): #loop over chunks to block average
        print('\tNow on chunk: '+str(i+1)+
              ' out of '+str(split)+'\n')
        for j in range(tn): #loop over timesteps in block
            for k in range(9): #skip comments
                fid.readline()
            for k in range(num): #get atoms
                tmp = fid.readline().strip().split()
                vels[j,k,0] = float(tmp[2]) #vx
                vels[j,k,1] = float(tmp[3]) #vy
                vels[j,k,2] = float(tmp[4]) #vz
        vels = np.fft.fft(vels,axis=0)*dt*dn #take fft along time axis
        velsFFT = velsFFT+vels 
        
velsFFT = velsFFT/float(split) #average across chunks    
velsFFT = mass/2.0*np.multiply(abs(velsFFT),abs(velsFFT))/(tn*dt*dn) #VPS
velsFFT = velsFFT.sum(axis=2) #dot product, vx*vx + vy*vy + vz*vz
velsFFT = velsFFT.mean(axis=1) #average across all atoms

del velsFile, vels, split, prints, steps, dt, dn, tn, mass, tmp
####

### PLOT RESULTS ###
velsFFT = ty.gsmooth(velsFFT,win,dom) #gaussian smoothing

plt.plot(thz,velsFFT,c='k',linewidth=2)
plt.xlabel('Frequency, THz')
plt.ylabel('VPS, arb. units')
plt.title('Vibrational Power Spectrum, Argon-like LJ Solid')
plt.axis([0,0.5,0,2.5e-17])
plt.show()
####

del a, dom, i, j, k, num, pi, win
    
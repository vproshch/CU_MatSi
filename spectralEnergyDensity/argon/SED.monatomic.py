#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 25 14:58:06 2018

@author: Ty Sterling

LICENSE AGREEMENT: Do whatever you want :)

i wrote all this in using python 2.7 installed with conda.
"""

### GLOBAL VARIABLES ###
import numpy as np #if it doesnt work, lemme know and ill tell you the versions
import matplotlib.pyplot as plt
import ty


velsFile = 'vels.dat'
posFile = 'data.argon'

steps = 50000 #run time
dt = 5e-15 #lammps time step
dn = 100 #print frequency
prints = steps/dn #times data is printed
split = 1 #times to split data for averaging
tn = prints/split #timesteps per chunk
win = 0.05 #gaussian smoothing window
pi = np.pi #tired of forgetting the 'np' part...

mass = 39.95 #argon
a = 5.376 #lattice constant of argon

prim = np.array([[0,0.5,0.5],
                 [0.5,0,0.5],
                 [0.5,0.5,0]])*a #primitive lattice vectors

specialk = np.array([[0,0,0], #gamma
                     [0.5,0,0.5], #X
                     [0.375,0.375,0.75], #K
                     [0,0,0], #gamma
                     [0.5,0.5,0.5]]) #L 
                     #special reciprocal lattice points
klabel = np.array(('G','X','K','G','L')) 
dk = 400 #k space mesh, number of points between speciak k points

kpoints, kdist = ty.makeKpoints(prim,specialk,dk) #get the input k space arrays
#from a funtion

### GET COORDS ###
num, types, masses, pos = ty.readData(posFile) #get the positions from the 
#lammps data file
####

## GET VELOCITIES AND CALCULATE SED ###
with open(velsFile, 'r') as fid: 
    nk = len(kpoints) #number of k points
    ids = np.zeros((num))
    
    sed = np.zeros((tn,nk)).astype(complex) #spectral energy density
    
    #the data is read in in chunks and the chunks, mathed upon until
    #its a smaller data structure, and the chunkns are all block averaged
    #together. Saves RAM space (trust me, i've got 16 GB's and I've maxed
    #it out doing transmission calculateions (20 GB of veclocity data!)) 
    #and it also 'ensemble averages' to produce better data                
    for i in range(split): #loop over chunks to block average
        print('\tNow on chunk: '+str(i+1)+
              ' out of '+str(split)+'\n')
        vels = np.zeros((tn,num,3))
        qdot = np.zeros((tn,nk))
        
        for j in range(tn): #loop over timesteps in block
            tmpVels = np.zeros((num,3))
            for k in range(9): #skip comments
                fid.readline()
            for k in range(num): #get atoms
                tmp = fid.readline().strip().split()
                vels[j,k,0] = float(tmp[2]) #vx
                vels[j,k,1] = float(tmp[3]) #vy
                vels[j,k,2] = float(tmp[4]) #vz
                
        for j in range(nk): #loop over all atoms
            kvec = kpoints[j,:] #k space vector
            tmp = np.zeros(tn)
            for k in range(num): #loop over x,y,z 
                rvec = pos[k,2:5] #real space vector
                tmp = (tmp+(vels[:,k,0]+vels[:,k,1]+vels[:,k,2])
                *np.exp(1j*(np.dot(kvec,rvec)))) 
            qdot[:,j] = np.multiply(abs(np.fft.fft(tmp[:])),
                abs(np.fft.fft(tmp[:]))) 
            #spectral energy density
        
        sed = sed+qdot*(mass/num)/(4*np.pi*dt*tn*dn)
        
sed = sed/split

### PLOT THE DISPERSION CURVE ###
plt.imshow(np.real(sed),cmap='plasma')
plt.tight_layout()
plt.show()
####

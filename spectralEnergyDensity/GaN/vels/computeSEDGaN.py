#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 25 14:58:06 2018

@author: Ty Sterling

ty.sterling@colorado.edu

LICENSE AGREEMENT: Do whatever you want :)

i wrote all this in using python 2.7 installed with conda.
"""

### GLOBAL VARIABLES ###
import numpy as np #if it doesnt work, lemme know and ill tell you the versions
import matplotlib.pyplot as plt
import ty

ty.tic()

outfile = 'sed.ty.GaN.dat'
velsfile = 'vels.dat'

steps = 50000 #run time
dt = 0.2e-15 #lammps time step
dn = 10 #print frequency
prints = steps/dn #times data is printed
split = 4 #times to split data for averaging
tn = prints/split #timesteps per chunk
win = 0.1 #gaussian smoothing window
pi = np.pi #tired of forgetting the 'np' part...

#om = np.arange(0,tn)*2*np.pi/(tn*dt*dn) #angular frequency
thz = np.arange(0,tn)/(tn*dt*dn)*1e-12 #frequency in THz

#### GET POSITIONS AND UNITCELL INDICIES ###
num, pos, masses, uc, a, c = ty.makeGaN(10,10,14) #get the positions from the 
#lammps data file
types = pos[:,1]
###

#### GET K POINTS ###
prim = np.array([[1,0,0],
                 [0,1,0],
                 [0,0,1]]).astype(float) #primitive lattice vectors

prim[0,:] = prim[0,:]*np.sqrt(3.0)*a #ax
prim[1,:] = prim[1,:]*a #ay
prim[2,:] = prim[2,:]*c #az

specialk = np.array([[0, 0, 0], #G
                     [0.5, 0, 0], #X
                     [0.5, 0.5, 0], #S
                     [0 ,0.5, 0], #Y
                     [0, 0, 0], #G
                     [0, 0, 0.5], #Z
                     [0.5, 0, 0.5], #U
                     [0.5, 0.5, 0.5], #R
                     [0, 0.5, 0.5], #T
                     [0, 0, 0.5]]) #Z 
                     #special reciprocal lattice points
klabel = np.array(('G','X','S','Y','G','Z','U','R','T','Z')) 
dk = 50 #k space mesh, number of points between speciak k points

kpoints, kdist = ty.makeKpoints(prim,specialk,dk) #get the input k space arrays
#from a funtion

### GET VELOCITIES AND CALCULATE SED ###
with open(velsfile, 'r') as fid: 
    nk = len(kpoints) #number of k points
    ids = np.zeros((num))
    nc = max(uc)+1
    nb = len(masses)
    ids = np.argwhere(types == 1) #atoms in this basis pos
    cellvec = pos[ids[:,0],:] #coords of fcc basis atom
    cellvec = cellvec[np.argsort(uc[ids][:,0]),:] #sort by unit cell
    
    sed = np.zeros((tn,nk)).astype(complex) #spectral energy density
    
    #the data is read in in chunks and the chunks are mathed upon until
    #its a smaller data structure. then the chunks are all block averaged
    #together. Saves RAM space \and it also 'ensemble averages' to 
    #produce better data      
          
    for i in range(split): #loop over chunks to block average
        print('\n\tNow on chunk: '+str(i+1)+
              ' out of '+str(split)+'\n')
        vels = np.zeros((tn,num,3))
        qdot = np.zeros((tn,nk))
        
        #read in vels for this block
        for j in range(tn): #loop over timesteps in block
            tmpVels = np.zeros((num,3))
            for k in range(9): #skip comments
                fid.readline()
            for k in range(num): #get atoms
                tmp = fid.readline().strip().split()
                vels[j,k,0] = float(tmp[2]) #vx
                vels[j,k,1] = float(tmp[3]) #vy
                vels[j,k,2] = float(tmp[4]) #vz
                
        #compute SED for this block        
        for j in range(nk): #loop over all k points
            kvec = kpoints[j,:] #k point vector
            if j%50 == 0:
                print('\t\tNow on k-point: '+str(j)+' out of '+str(nk))
            tmp = np.zeros((tn,1)).astype(complex) #sed for this k point
            for k in range(nc-1): #loop over unit cells
                rvec = cellvec[k,2:5] #position of unit cell
                ids = np.argwhere(uc==k) 
                for l in range(nb): #loop over basis atoms
                    mass = masses[l]
                    vx = vels[:,ids[l,0],0] #time series for particular atom
                    vy = vels[:,ids[l,0],1] # '' ''
                    vz = vels[:,ids[l,0],2] # '' ''
                    
                    tmp[:,0] = tmp[:,0]+(vx+vy+vz)*np.exp(1j*np.dot(kvec,rvec))
                    #space fft of velocity data
                    qdot[:,j] = np.multiply(abs(np.fft.fft(tmp[:,0])),
                                      abs(np.fft.fft(tmp[:,0])))*mass 
                    #KE of normal coordinate (square of time-FFT)
                    
        sed = sed+qdot/(4*np.pi*steps/split*dt*nc) #a buncha constants

sed = sed/split #average across splits

### WRITE TO A FILE ###
ty.writeSED(outfile,thz,kpoints,sed)

ty.toc()

### PLOT THE DISPERSION CURVE ###
plt.imshow(np.real(sed),cmap='plasma',aspect='auto')
plt.tight_layout()
plt.show()
###


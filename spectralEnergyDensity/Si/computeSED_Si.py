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

outfile = 'sed.ty.SHORT.dat'
velsfile = 'vels.dat'

steps = 500000 #run time
dt = 0.5e-15 #lammps time step
dn = 40 #print frequency
prints = steps/dn #times data is printed
split = 5 #times to split data for averaging
tn = prints/split #timesteps per chunk
win = 0.05 #gaussian smoothing window
pi = np.pi #tired of forgetting the 'np' part...

nx = 16; ny = 16; nz = 16 #size of structure in x, y, and z (# of unit cells)
dk = 40 #k space mesh, number of points between speciak k points

#om = np.arange(0,tn)*2*np.pi/(tn*dt*dn) #angular frequency
thz = np.arange(0,tn)/(tn*dt*dn)*1e-12 #frequency in THz

### GET POSITIONS AND UNITCELL INDICIES ###
num, pos, masses, uc, a = ty.makeFCCdiamond(nx,ny,nz) #get the positions from 
#the lammps data file
types = pos[:,1]
###

### GET K POINTS ###
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

kpoints, kdist = ty.makeKpoints(prim,specialk,dk) #get the k point vectors
#from a funtion

## GET VELOCITIES AND CALCULATE SED ###
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
    #together. Saves RAM space (trust me, i've got 16 GB's and I've maxed
    #it out doing transmission calculateions (20 GB of veclocity data!)) 
    #and it also 'ensemble averages' to produce better data      
          
    for i in range(2):#split): #loop over chunks to block average
        print('\n\tNow on chunk: '+str(i+1)+
              ' out of '+str(split)+'\n')
        vels = np.zeros((tn,num,3))
        qdot = np.zeros((tn,nk))
        t = steps/split*dt
        
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
            for m in range(1): #sum over x, y, z
                for k in range(nb): #loop over basis atoms
                    mass = masses[k] #mass of this basis atom
                    ids = np.argwhere(types == k+1)
                    ids = ids[np.argsort(uc[ids][:,0])] #sort by unit cell
                    for l in range(nc): #sum over all unit cells
                        atom = ids[l] #particular atom
                        rvec = cellvec[l,2:5]
                        vt = vels[:,atom,m] #time series for particular atom
                        tmp[:] = tmp[:]+(vt)*np.exp(1j*np.dot(kvec,rvec))
                        #space fft of velocity data
                        qdot[:,j] = np.multiply(abs(np.fft.fft(tmp[:,0])),
                                          abs(np.fft.fft(tmp[:,0])))*mass 
                        #KE of normal coordinate (square of time-FFT)
                    
        sed = sed+qdot#/(4*np.pi*steps/split*dt*nc) #a buncha constants

#sed = sed/split #average across splits

### WRITE TO A FILE ###
ty.writeSED(outfile,thz,kpoints,sed)

### PLOT THE DISPERSION CURVE ###
#plt.imshow(np.real(sed),cmap='jet',aspect='auto')
#plt.tight_layout()
#plt.show()
####

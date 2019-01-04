#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Calculate phonon band structure by diagonalizing the dynamical matrix
Force constants are calculated in LAMMPS
"""
import numpy as np
import matplotlib.pyplot as plt
import ty

### INPUT DATA ###
posfile = 'data.si'
forcefile = 'Fij.dat'

a = 5.431 #silicon lattice constant
basis = np.array((274,306)) #basis atoms, pick 2 in the middle somewhere
prim = np.array([[0,0.5,0.5],
                 [0.5,0,0.5],
                 [0.5,0.5,0]])*a #primitive lattice vectors
diamond = a/4.0 #tranlation vector to 2nd basis atom

specialk = np.array([[0,0,0], #gamma
                     [0.5,0,0.5], #X
                     [0.375,0.375,0.75], #K
                     [0,0,0], #gamma
                     [0.5,0.5,0.5]]) #L 
                     #special reciprocal lattice points
klabel = np.array(('G','X','K','G','L')) 
dk = 100 #k space mesh, number of points between speciak k points
###

### CALCULATE K-SPACE POINTS ###
kpoints, kdist = ty.makeKpoints(prim,specialk,dk)
del prim, a
####

### READ IN COORDINATES OF ATOMS ###
num, types, masses, pos = ty.readData(posfile) #number of atoms, number of 
#atom types, mass of each type, and positions vector [:,(ID, TYPE, X, Y, Z)]
masses = np.array([28.0855,28.0855])
del posfile
####

### FIND NEAREST NEIGHBOR LIST FOR BASIS ATOMS ###
nl, dl, nn = ty.findNN(pos)
####

### READ IN FORCE CONSTANT FROM LAMMPS ###
phi, rvec = ty.readFij(forcefile,num,pos,basis,nl,nn,neighbor=2,tol=1e-5) 
#force constant matrix for each basis and neighbor atom, direction vector
#between each basis and atom and neighbor atom
del forcefile 
####

### COMPUTE THE DYNAMICAL MATRIX AND DIAGONALIZE ###
om = ty.dynamicalMat(phi,rvec,pos,basis,masses,kpoints,nl)
####
    
### PLOT DISPERSION ###
ommax = om.max(axis=0).max() #used for plotting

plt.plot(np.ones((100,1))*kdist[0], #plot solid veritcal bars
         np.linspace(0,ommax*1.1,100),'k-',linewidth=2)
 
for i in range(len(klabel)-1):
    plt.plot(np.linspace(kdist[i],kdist[i+1],dk),om[i*dk:(i+1)*dk,:],
             'm.',markersize=2)
    
    plt.plot(np.ones((100,1))*kdist[i+1],
         np.linspace(0,ommax*1.1,100),'k-',linewidth=2)
    
plt.xticks((kdist[0],kdist[1],kdist[2],kdist[3],kdist[4]),klabel)
plt.ylabel('omega, THz')

plt.show()





        
        
            

            
            
            

    



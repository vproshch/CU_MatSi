#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Calculate and Plot dispersion curves for bulk crystalline Silicon
Based on MATLADY codes by brucefan1983 on github.

This code is not generic; it only works for FCC silicon.
Process is:
    1: Create silicon FCC unit cell of 8 atoms
        Array to create bulk strucutre.
        Select some atoms representing an fcc site and a diamond site to
        form the basis; pick these atom in the center of the 4x4x4 cubic 
        structure so that they are far from the boundaries.
    2: From primitive translation vectors and knowledge of silicon reciprocal
        lattice, create k-space array between special points in k space
        i.e. from gamma to K, K to X, X to gamma, and gamma to L
    3: Find each atoms nearest neighbors. The only relevant atoms are the
        basis atoms, thier nearest neighbors, and thier neighbors nearest
        neighbors. We need to know these to compute the forces constants.
    4: Calculate force constants by finite different derivatives, i.e.
        move some atoms a small amount and calculate the change in lattice
        energy. Assuming small displacements, the 2nd order term in the 
        Taylor expansion of the lattice is the force constant, 
        i.e. d^2V/dui_a/duj_b. Each force constant is a 3x3 matrix where
        element i,j = 1:3 is the force constant component in the x,y,z 
        direction.
    5: The force constants are calculated from the tersoff potential.
        For each force contant between a basis atom and it's nearest neighbors,
        the lattice energy is iteratively calculated with the basis and 
        neighbor atoms displaced slightly. The difference in lattice energies,
        d^2V, is calculated. The ratio of of the difference in energy to the 
        displacement is the force constant.   
    6: The dynamical matrix is constructed from the force constant matricies.
        The dynamical matrix is formed by creating the hessian force constant
        matrix about each basis atom and its nearest neigbors. e.g. for Si.
        the force constant matrix is 6x6 corresponding to dV^2/dui_a/duj_b
        where each i,j are the force contants for basis atoms 1 in x,y,z '
        and 2 in x,y,z with neighbor atoms 1 in x,y,z and 2 in x,y,z. 
        The dynamical matrix is formed by Fourier tranforming the force 
        constant matricies. The frequencies and displacement vectors are 
        found by diagonalizing the dynamical matrix.
"""

import numpy as np
import matplotlib.pyplot as plt
import copy as cp
pi = np.pi

### DEFINE STRUCTURE ###
a = 5.431 #Si lattice constant
scale = a/4.0 #multiply final positions by this
nT = np.array([4,4,4]) #times to translate cell in each direction
edges = nT*a #edge lengths of structure

unit = np.array(([0,0,0,0], #FCC-diamond unit cell (basis-index,x,y,z)
                  [0,0,2,2], #Index 0 are FCC sites, 1 are diamond sites
                  [0,2,0,2],
                  [0,2,2,0],
                  [1,1,1,1],
                  [1,1,3,3],
                  [1,3,1,3],
                  [1,3,3,1])).astype(float)

basis = np.array([336,340]) #basis atoms to find force constants
masses = np.array([28.0855,28.0855]) #atomic mass of Si
####

### DEFINE POTENTIAL ###
tersoff = np.array((1.8308e3,471.18,2.4799,1.7322,1.1e-6,0.78734,1.0039e5,
                    16.217,-0.59825,2.7,3)) #SiCGe.tersoff parameters
####

### CREATE STRUCTURE ###
tmp = cp.deepcopy(unit)
pos = cp.deepcopy(unit) #basis-index and coordinates of atoms
for i in range(nT[0]-1): #translate in x-direction 
    tmp[:,1] = tmp[:,1]+4
    pos = np.append(pos,tmp,axis=0)

tmp = cp.deepcopy(pos)
for i in range(nT[1]-1): #translate in y-direction
    tmp[:,2] = tmp[:,2]+4
    pos = np.append(pos,tmp,axis=0)

tmp = cp.deepcopy(pos)
for i in range(nT[2]-1): #translate in z-direction
    tmp[:,3] = tmp[:,3]+4
    pos = np.append(pos,tmp,axis=0)
    
num = len(pos[:,0])
ids = np.zeros((num,1)) 
ids[:,0] = np.arange(0,num,1) #will be used later

pos = np.append(ids,pos,axis=1)
unit = pos[0:len(unit[:,0])+1,:]
unit = unit[:,2:5]*scale #rescale coordinates
pos[:,2:5] = pos[:,2:5]*scale

del tmp, i, nT, ids
####

### GENERATE K-POINTS ###
prim = np.array(([0,1,1],[1,0,1],[1,1,0]))*a/2.0 #primitive lattice vectors
specialK = np.array(([0,0,0], #gamma
                     [1/2.0,0,1/2.0], #X
                     [3/8.0,3/8.0,3/4.0], #K
                     [0,0,0], #gamma
                     [1/2.0,1/2.0,1/2.0])) #L
klabel = np.array(('G','X','K','G','L')) #labels for special K-points
dK = 100 #mesh size between special K-points

recip = np.zeros((3,3)) #reciprocal lattice vectors
vol = np.dot(prim[:,0],np.cross(prim[:,1],prim[:,2])) #primitive cell volume

recip[:,0] = 2*pi*np.cross(prim[:,1],prim[:,2])/vol #reciprocal lattice vector
recip[:,1] = 2*pi*np.cross(prim[:,2],prim[:,0])/vol #see, e.g. Ashcroft-Mermin
recip[:,2] = 2*pi*np.cross(prim[:,0],prim[:,1])/vol

for i in range(len(specialK)): #reciprocal space coordinates of K-points 
    specialK[i,:] = np.dot(recip,specialK[i,:]) 
    
nK = (len(specialK)-1)*dK #len of K-points array
kPoints = np.zeros((nK,3))

for i in range(len(specialK)-1): #populate the array between special K-points
    for j in range(3):
        kPoints[i*dK:(i+1)*dK,j] = np.linspace(specialK[i,j],specialK[i+1,j],
                num=len(kPoints[i*dK:(i+1)*dK,j]))

kDist = np.zeros((len(specialK)))
kDist[0] = 0
for i in range(len(specialK)-1):
    kDist[i+1] = np.linalg.norm(specialK[i+1,:]-specialK[i,:]) 
kDist = np.cumsum(kDist) #total distance 'K-vector' travels through K-space

del prim, specialK, nK, vol, recip, i, j 
####

### FIND NEAREST NEIGHBOR LISTS WITH PBC ###
nn1 = np.zeros((num,1)).astype(int) #number of 1st nearest neighbors
nn2 = np.zeros((num,1)).astype(int) #number of 2nd nearest neighbors
nl1 = np.zeros((num,num-1)).astype(int) #1st neigbor list (excluding self)
nl2 = np.zeros((num,num)).astype(int) #1st+2nd neighbor list (including self)
dl1 = np.zeros((num,num-1)) #1st nn distances
dl2 = np.zeros((num,num)) #2nd nn distances
 
for i in range(num):
    tmp = np.zeros((num,2))
    tmp[:,0] = pos[:,0] #atom ID's
    for j in range(num):
        tmp[j,1] = np.round(np.sqrt((pos[i,2]-pos[j,2])**2+
           (pos[i,3]-pos[j,3])**2+(pos[i,4]-pos[j,4])**2),decimals=12)
    tmp = tmp[np.argsort(tmp[:,1]),:] #sort list to find nn's
    nn1[i] = sum(tmp[:,1] == tmp[1,1]) #skip self
    nn2[i] = sum(tmp[:,1] == tmp[int(nn1[i,0])+1,1])+nn1[i,0]+1
    nl2[i,:] = tmp[:,0] #2nd nn's
    dl2[i,:] = tmp[:,1] #2nd nn distances
    tmp = np.delete(tmp,0,axis=0) #remove self from list
    nl1[i,:] = tmp[:,0] #1st nn's
    dl1[i,:] = tmp[:,1] #1st nn distances
    
nl1 = nl1[:,0:4] #all silicon sites have 4 1st nn's
nl2 = nl2[:,0:17] #all silicon sites have 12 2nd nn's (1st nn+2nd nn+self)
dl1 = dl1[:,0:4] #1st nn distance
dl2 = dl2[:,0:17] #2nd nn distance

del i, j, tmp
####

### FIND FORCE CONSTANTS ###
###Take finite difference derivative of lattice energy by shift atoms

A = tersoff[0]; B = tersoff[1]; lam = tersoff[2]; mu = tersoff[3]; 
beta = tersoff[4]; eta = tersoff[5]; c = tersoff[6]; d = tersoff[7];
h = tersoff[8]; r1 = tersoff[9]; r2 = tersoff[10]; c2 = c**2; d2 = d**2
opcod = 1+c2/d2; piFactor = pi/(r2-r1); mhon = -0.5/eta
#parameters for force calculation

phi = [0]*len(basis) #force constant for each basis atom
dR = [0]*len(basis) #direction vector between basis atom and nn's

for i in range(len(basis)): #populate force constant matricies
    atom = basis[i] #translational invariance, only need to check basis atoms
    phi[i] = [0]*int(nn2[atom]) #nn2 no. of force constant matricies
    dR[i] = [0]*int(nn2[atom]) #direction vectors to nn's
    for j in range(nn2[atom]): #loop over nn's
        hess = np.zeros((3,3)) #hessian force constant matrix
        du = 0.0025 #force constant = d^2V/du^2
        for k in range(3): #x, y, z 
            for l in range(3): #x, y, z for mixed derivatives
                rpp = cp.deepcopy(pos) #shifted positions in +x,y,z; +x,y,z
                rmm = cp.deepcopy(pos) #-x,y,z; -x,y,z
                rpm = cp.deepcopy(pos) #+x,y,z; -x,y,z
                rmp = cp.deepcopy(pos) #-x,y,z; +x,y,z
                
                rpp[atom,k+2] = rpp[atom,k+2]+du #shift basis atom
                rpp[nl2[atom,j],l+2] = rpp[nl2[atom,j],l+2]+du #shift nn
                rmm[atom,k+2] = rmm[atom,k+2]-du #-x,y,z
                rmm[nl2[atom,j],l+2] = rmm[nl2[atom,j],l+2]-du #-x,y,z
                rpm[atom,k+2] = rpm[atom,k+2]+du #+x,y,z
                rpm[nl2[atom,j],l+2] = rpm[nl2[atom,j],l+2]-du #-x,y,z
                rmp[atom,k+2] = rmp[atom,k+2]-du #-x,y,z
                rmp[nl2[atom,j],l+2] = rmp[nl2[atom,j],l+2]+du #+x,y,z

                data = [[atom]+[nl2[atom,j]]+
                        nl1[atom,0:int(nn1[atom])].tolist()+
                        nl1[nl2[atom,j],0:int(nn1[nl2[atom,j]])].tolist()]
                #basis site, basis site's nn, all of basis site's
                #1st nn's, all of basis site's nn's 1st nn's
                #lattice sites used to calculate lattice energy
                data = np.array(data)
                data = np.unique(data) #unique list of lattice sites with PBC
                
                ### RPP 
                tmp = rpp #temporary positions variable
                energy = 0
                for m in range(len(data)): #sum over lattice sites
                    site1 = int(data[m]) #central lattice site
                    for n in range(nn1[site1]): #sum over sites nn's
                        dist12 = dl1[site1,n] #distance between site and nn
                        site2 = nl1[site1,n]
                        r12 = tmp[site2,2:5]-tmp[site1,2:5]
                     
                        fR = A*np.exp(-lam*dist12)
                        fA = B*np.exp(-mu*dist12)
                        if dist12 < r1: #cut off distance
                            fC = 1.0
                        else:
                            fC = np.cos(piFactor*(dist12-r1))/2.0+1/2.0
                        zeta = 0.0
                        for p in range(nn1[site1]):
                            site3 = nl1[site1,p]
                            if site3 == site2:
                                continue   
                            dist13 = dl1[site1,p]
                            r13 = tmp[site3,2:5]-tmp[site1,2:5]
                            
                            cos123 = np.dot(r12,r13)/(dist12*dist13)
                            if dist13 < r1: #cut off distance
                                fC13 = 1.0
                            else:
                                fC13 = np.cos(piFactor*(dist13-r1))/2.0+1/2.0
                            g = opcod-c2/(d2+(cos123-h)*(cos123-h))
                            zeta = zeta+fC13*g
                        b12 = (1+(beta*zeta)**eta)**mhon
                        energy = energy+fC/2.0*(fR-b12*fA)
                epp = energy
                ###
                
                ### RMM
                tmp = rmm #temporary positions variable
                energy = 0
                for m in range(len(data)): #sum over lattice sites
                    site1 = int(data[m]) #central lattice site
                    for n in range(nn1[site1]): #sum over sites nn's
                        dist12 = dl1[site1,n] #distance between site and nn
                        site2 = nl1[site1,n]
                        r12 = tmp[site2,2:5]-tmp[site1,2:5]
                     
                        fR = A*np.exp(-lam*dist12)
                        fA = B*np.exp(-mu*dist12)
                        if dist12 < r1: #cut off distance
                            fC = 1.0
                        else:
                            fC = np.cos(piFactor*(dist12-r1))/2.0+1/2.0
                        zeta = 0.0
                        for p in range(nn1[site1]):
                            site3 = nl1[site1,p]
                            if site3 == site2:
                                continue   
                            dist13 = dl1[site1,p]
                            r13 = tmp[site3,2:5]-tmp[site1,2:5]
                            
                            cos123 = np.dot(r12,r13)/(dist12*dist13)
                            if dist13 < r1: #cut off distance
                                fC13 = 1.0
                            else:
                                fC13 = np.cos(piFactor*(dist13-r1))/2.0+1/2.0
                            g = opcod-c2/(d2+(cos123-h)*(cos123-h))
                            zeta = zeta+fC13*g
                        b12 = (1+(beta*zeta)**eta)**mhon
                        energy = energy+fC/2.0*(fR-b12*fA)
                emm = energy
                ###
                
                ### RPM
                tmp = rpm #temporary positions variable
                energy = 0
                for m in range(len(data)): #sum over lattice sites
                    site1 = int(data[m]) #central lattice site
                    for n in range(nn1[site1]): #sum over sites nn's
                        dist12 = dl1[site1,n] #distance between site and nn
                        site2 = nl1[site1,n]
                        r12 = tmp[site2,2:5]-tmp[site1,2:5]
                     
                        fR = A*np.exp(-lam*dist12)
                        fA = B*np.exp(-mu*dist12)
                        if dist12 < r1: #cut off distance
                            fC = 1.0
                        else:
                            fC = np.cos(piFactor*(dist12-r1))/2.0+1/2.0
                        zeta = 0.0
                        for p in range(nn1[site1]):
                            site3 = nl1[site1,p]
                            if site3 == site2:
                                continue   
                            dist13 = dl1[site1,p]
                            r13 = tmp[site3,2:5]-tmp[site1,2:5]
                            
                            cos123 = np.dot(r12,r13)/(dist12*dist13)
                            if dist13 < r1: #cut off distance
                                fC13 = 1.0
                            else:
                                fC13 = np.cos(piFactor*(dist13-r1))/2.0+1/2.0
                            g = opcod-c2/(d2+(cos123-h)*(cos123-h))
                            zeta = zeta+fC13*g
                        b12 = (1+(beta*zeta)**eta)**mhon
                        energy = energy+fC/2.0*(fR-b12*fA)
                epm = energy
                ###
                
                ### RMP
                tmp = rmp #temporary positions variable
                energy = 0
                for m in range(len(data)): #sum over lattice sites
                    site1 = int(data[m]) #central lattice site
                    for n in range(nn1[site1]): #sum over sites nn's
                        dist12 = dl1[site1,n] #distance between site and nn
                        site2 = nl1[site1,n]
                        r12 = tmp[site2,2:5]-tmp[site1,2:5]
                     
                        fR = A*np.exp(-lam*dist12)
                        fA = B*np.exp(-mu*dist12)
                        if dist12 < r1: #cut off distance
                            fC = 1.0
                        else:
                            fC = np.cos(piFactor*(dist12-r1))/2.0+1/2.0
                        zeta = 0.0
                        for p in range(nn1[site1]):
                            site3 = nl1[site1,p]
                            if site3 == site2:
                                continue   
                            dist13 = dl1[site1,p]
                            r13 = tmp[site3,2:5]-tmp[site1,2:5]
                            
                            cos123 = np.dot(r12,r13)/(dist12*dist13)
                            if dist13 < r1: #cut off distance
                                fC13 = 1.0
                            else:
                                fC13 = np.cos(piFactor*(dist13-r1))/2.0+1/2.0
                            g = opcod-c2/(d2+(cos123-h)*(cos123-h))
                            zeta = zeta+fC13*g
                        b12 = (1+(beta*zeta)**eta)**mhon
                        energy = energy+fC/2.0*(fR-b12*fA)
                emp = energy
                ###

                hess[k,l] = (epp+emm-epm-emp)/du/du/4.0 #force constant mat.
        phi[i][j] = hess #write fc matrix to array
        dR[i][j] = pos[nl2[atom,j],2:5]-pos[atom,2:5]
        
del A, B, b12, atom, beta, c, c2, d, d2, cos123, data, dist12, dist13
del du, emm, epp, emp, epm, fA, fC, fC13, fR, g, h, hess, i, j, k, l
del lam, energy, eta, m, mhon, mu, n, opcod, p, piFactor, r1, r12, r13, r2
del rmm, rmp, rpm, rpp, site1, site2, site3, tersoff, tmp, zeta

#### CALCULATE DIAGONAL MATRIX, 
#nb = len(basis[:]) #number of basis atoms
#nk = len(kPoints[:,0]) #number of k-points
#omega = np.zeros((nk,nb*3)) #frequencies, 3 solutions per basis atom: 3 accoutic 
#                         #modes, 3*(nb-1) optic modes.
#                         
#for i in range(nk):
#    kvec = kPoints[i,:] #solve at a particualr k point
#    dyn = np.zeros((nb*3,nb*3)).astype(complex) #dynamical matrix
#    
#    for j in range(nb):
#        atom1 = basis[j] #ID of basis atom
#        type1 = int(pos[atom1,1]) #type of basis atom
#        mass1 = masses[type1] #mass of basis atom
#        index1 = np.arange(0,3)+type1*3 #dx,dy,dz in dyn for basis atom
#        for k in range(nn2[atom1]):
#            atom2 = nl2[atom1,k] #ID of neighbor atom
#            type2 = int(pos[atom2,1]) #type of neighbor atom
#            mass2 = masses[type2] 
#            index2 = np.arange(0,3)+type2*3 #dx,dy,dz in dyn for neighbor atom
#            dyn[index1[0]:index1[-1]+1,index2[0]:index2[-1]+1] = (
#                    dyn[index1[0]:index1[-1]+1,index2[0]:index2[-1]+1]+
#                    phi[j][k]*np.exp(1j*np.dot(kvec,dR[j][k]))
#                    /np.sqrt(mass1*mass2)) #space FFT of force constant matrix
#                    #dived by the sqrt of the product of all basis masses
#    eigVal, eigVec = np.linalg.eig(dyn) #eigen values are frequency, vectors
#    #are displacement vectors
#    omega[i,:] = np.sqrt(np.real(eigVal))*1000/2.0/pi/10.18
#
#del atom1, type1, mass1, index1, atom2, type2, mass2
#del index2, i, j, k, kvec, a, dl1, dl2, edges, masses, nb, nk    
#####
#    
#### PLOT THE DISPERSION CURVE ###
#omMax = omega.max(axis=0).max() #used for plotting
#plt.plot(np.ones((100,1))*kDist[0], #plot solid veritcal bars
#         np.linspace(0,omMax*1.1,100),'k-',linewidth=2) 
#for i in range(len(klabel)-1):
#    plt.plot(np.linspace(kDist[i],kDist[i+1],dK),omega[i*dK:(i+1)*dK,:],
#             'm.',markersize=2)
#    plt.plot(np.ones((100,1))*kDist[i+1],
#         np.linspace(0,omMax*1.1,100),'k-',linewidth=2)
#plt.xticks((kDist[0],kDist[1],kDist[2],kDist[3],kDist[4]),klabel)
#plt.ylabel('omega, THz')
#####
            
        
        
        
    
    


#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 20 17:12:28 2018

@author: ty

For fcc diamond lattice, unit cell coords are:
[[0,0,0],[0,0,a],[0,a/2,a/2],[0,a,0],[0,a,a],[a/2,0,a/2],[a/2,a/2,0],[a/2,a/2,a],[a/2,a,a/2],[a,0,0],
[a,0,a],[a,a/2,a/2],[a,a,0],[a,a,a],[a/4,a/4,a/4],[3*a/4,3*a/4,a/4],[a/4,3*a/4,3*a/4],[3*a/4,a/4,3*a/4]]

"""   
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import copy as cp
import sys
import random
import math




################    GET INPUT FROM USER    #####################
print ('\n\n\nThis program initializes FCC diamond structures, e.g. C-diamond, Si, Ge, and Si/Ge heterostrucutres.\n'
       'For Si and Si/Ge, the script initializes the whole structure using the Si lattice constant, a=5.431 A.\n'
       'For C, a=3.57 A. For pure Ge, a=5.658 A. For Si/Ge, the "crossplane" direction is the x-direction\n'
       '(i.e. nx = (# of Si unit cells) + (# of Ge unit cells)) with Si on the "left."\n\n'
       'Additionally, defects can be included in the structure, e.g. subsitutionals, self-interstitials,\n'
       'impurity interstitials, vacancies, vacancy pairs, and vacancy-interstitial pairs (Frenkel Pairs).\n\n' 
       'Defects are placed randomly and uniformly and can be chosen to populate the entire structure or to populate\n'
       'only a range specified spanning from left to right in the x-directions. This is useful for placing defects\n'
       'near an interface or for secluding defects from heat-baths on each end in NEMD.\n'  
       '------------------------------------------------------------------------------------------------\n\n')

###get species type from user
species = raw_input('Enter a species type for the bulk strcuture: (USAGE: C, Si, Ge, or Si/Ge)\n>')

###set lattice constant
if species == 'C':
    b = 3.57
elif species == 'Si':
    b = 5.431
elif species == 'Si/Ge':
    b = 5.431
else:
    b = 5.658
    
###check validity of input
if (species == 'C') or (species == 'Si') or (species == 'Ge') or (species == 'Si/Ge'):
    print ('\n\tYou entered: ' + species + '\n\tLattice constant is: ' + str(b) + ' A\n')
else:
    sys.exit('USAGE ERROR: Enter either C, Si, Si/Ge, or Ge')

###get size of material
if species == 'Si/Ge':
    n_si = int(raw_input('\nEnter the number of Si unit cells: (USAGE: an integer greater than 0)\n>'))
    n_ge = int(raw_input('\nEnter the number of Ge unit cells: (USAGE: an integer greater than 0)\n>'))
    if (n_si <= 0) or (n_ge <=0):
        sys.exit('USAGE ERROR: the number of unit cells must be greater than 0')
    nx = n_si + n_ge
    ny = int(raw_input('\nEnter number of unit cells in y-direction: (USAGE: an interger greater than 0)\n>'))
    if (ny <= 0):
        sys.exit('USAGE ERROR: the number of unit cells must be greater than 0')
    nz = int(raw_input('\nEnter number of unit cells in z-direction: (USAGE: an interger greater than 0)\n>'))
    if (nz <= 0):
        sys.exit('USAGE ERROR: the number of unit cells must be greater than 0')
else:
    nx = int(raw_input('\nEnter number of unit cells in x-direction: (USAGE: an interger greater than 0)\n>'))
    if (nx <= 0):
        sys.exit('USAGE ERROR: the number of unit cells must be greater than 0')
    ny = int(raw_input('\nEnter number of unit cells in y-direction: (USAGE: an interger greater than 0)\n>'))
    if (ny <= 0):
        sys.exit('USAGE ERROR: the number of unit cells must be greater than 0')
    nz = int(raw_input('\nEnter number of unit cells in z-direction: (USAGE: an interger greater than 0)\n>'))
    if (nz <= 0):
        sys.exit('USAGE ERROR: the number of unit cells must be greater than 0')
    
###display system size info    
print ('\n\tnx = ' + str(nx) + '\n\tny = ' + str(ny) + '\n\tnz = ' + str(nz) + '\n\n\tTotal Number of atoms = ' + 
           str(nx*ny*nz*8))

###get defected region and check input
defect = raw_input('\nAre there any defects? (USAGE: yes or no)\n>')
if (defect != 'yes') and (defect != 'no'):
    sys.exit('\nUSAGE ERROR: enter either yes or no')
if defect == 'yes':
    left_block = raw_input('\n\n------------------------------------------------------------------------------------\n'
                           'Select a range to include defects or type "all" to include defects in the entire structure. '
                           'If a range is specified,\ndefects are placed in and between the "left most" and "right most" '
                           'specified unit cells. For Si/Ge structures, input\nfor each side of the interface will '
                           'be collected later.\n\nLEFT BLOCK: (USAGE: "all" or an interger from 1 to ' + str(nx-1) + ')\n>')
    
    if left_block != 'all':
        all_flag = False
        try:
            left_block = int(left_block)
        except ValueError:
            sys.exit('USAGE ERROR: enter "all" or an integer from 1 to ' + str(nx-1) + '.')
        r = range(1,nx)
        if left_block not in r:
            sys.exit('USAGE ERROR: enter "all" or an integer from 1 to ' + str(nx-1) + '.')
        right_block = int(raw_input('\nRIGHT BLOCK: (USAGE: an integer from 2 to ' + str(nx) + ')\n>'))
        r = range(2,(nx+1))
        if right_block not in r:
            sys.exit('USAGE ERROR: please enter an integer from 2 to ' + str(nx) + '.')
        del r
    else:
        all_flag = True
        left_block = int(1)
        right_block = int(nx)
###########################################################





#############   MAKE STRUCTURE    ########################
###unit cell
a = 4

basis = np.array([[0,0,0],[0,a/2,a/2],[a/2,0,a/2],[a/2,a/2,0],[a/4,a/4,a/4],[3*a/4,3*a/4,a/4],
         [a/4,3*a/4,3*a/4],[3*a/4,a/4,3*a/4]])
intbasis = np.array([[a/4,a/4,3*a/4],[a/4,3*a/4,a/4],[3*a/4,a/4,a/4],[3*a/4,3*a/4,3*a/4]])    


###replicate in x direction
data = cp.deepcopy(basis)
intdata = cp.deepcopy(intbasis)
for i in range(1,nx):
    tmp = cp.deepcopy(basis)
    for j in range(len(tmp)):
        tmp[j][0] = tmp[j][0] + i*a
    data = np.append(data,tmp,axis=0)

new_basis = cp.deepcopy(data)
del basis

for i in range(1,nx):
    tmp = cp.deepcopy(intbasis)
    for j in range(len(tmp)):
        tmp[j][0] = tmp[j][0] + i*a
    intdata = np.append(intdata,tmp,axis=0)

intbasis = cp.deepcopy(intdata)

###replicate in y direction
for k in range(1,ny):
    tmp = cp.deepcopy(new_basis)
    for l in range(len(tmp)):
        tmp[l][1] = tmp[l][1] + k*a
    data = np.append(data,tmp,axis=0)
    
new_basis = cp.deepcopy(data)

for k in range(1,ny):
    tmp = cp.deepcopy(intbasis)
    for j in range(len(tmp)):
        tmp[j][1] = tmp[j][1] + k*a
    intdata = np.append(intdata,tmp,axis=0)

intbasis = cp.deepcopy(intdata)

###replicate in z direction
for m in range(1,nz):
    tmp = cp.deepcopy(new_basis)
    for n in range(len(tmp)):
        tmp[n][2] = tmp[n][2] + m*a
    data = np.append(data,tmp,axis=0)
    
new_basis = cp.deepcopy(data)    

for m in range(1,nz):
    tmp = cp.deepcopy(intbasis)
    for j in range(len(tmp)):
        tmp[j][2] = tmp[j][2] + m*a
    intdata = np.append(intdata,tmp,axis=0)

intbasis = cp.deepcopy(intdata)

###remove repeats
data = np.unique(new_basis,axis=0)
data = data[np.lexsort((data[:,2],data[:,1],data[:,0]))]
intdata = np.unique(intbasis,axis=0)
intdata = intdata[np.lexsort((intdata[:,2],intdata[:,1],intdata[:,0]))]
del new_basis, intbasis, tmp

ids = np.zeros((len(data),1))
types = np.zeros((len(data),1))

tmp = np.zeros((len(intdata),2))
intdata = np.append(tmp,intdata,axis=1)

for i in range(len(intdata)):
    intdata[i,0] = i+1
    intdata[i,1] = 1

###write types and ID's
if (species != 'Si/Ge'):
    for i in range(len(data)):
        if (species == 'C'):
            types[i] = 1
        elif (species == 'Si'):
            types[i] = 2
        else:
            types[i] = 3
        ids[i] = i+1
else:
    q = n_si*ny*nz*8
    r = n_si*ny*nz*4
    for j in range(q):
        types[j] = 2
        ids[j] = j+1
        
    for j in range(r):
        intdata[j,1] = 2
    for k in range(r,len(intdata)):
        intdata[k,1] = 3
        
    for k in range(q,len(data)):
        types[k] = 3
        ids[k] = k+1
        
data = np.append(types, data, axis=1)
data = np.append(ids, data, axis=1)
size = len(data)
del types, ids
####################################################





############   INCLUDE DEFECTS    ##################
###find out which defect types to include
if str(defect) != 'no':
    region = data[((left_block-1)*ny*nz*8):(right_block*ny*nz*8)][:]
    data = np.delete(data, range(((left_block-1)*ny*nz*8),(right_block*ny*nz*8)), axis=0)
    atoms = len(region)
    intdata = intdata[((left_block-1)*ny*nz*4):(right_block*ny*nz*4)][:]
#    intdata = np.delete(intdata, range(((left_block-1)*ny*nz*4),(right_block*ny*nz*4)), axis=0)
    print('\n\n\tThere are ' + str(len(region)) + ' Atoms in the defected region.\n')
    split_region = 'no'
    if (species == 'Si/Ge'):# and (all_flag != True):
        split_region = raw_input('\n-----------------------------------------------------------'
                                 '-----------------------------------------------------------\n'
                                 'This is a Si/Ge heterostructure with ' + str(n_si) + ' Si unit cell layers '
                                 'on the left and ' + str(n_ge) + ' Ge unit cell layers on the right.\n\n'
                                 'Are there different defects on each side of the interface? '
                                 '(USAGE: yes for "interfacial mixing", no for uniform defects)\n>')
        if (split_region != 'yes') and (split_region != 'no'):
            sys.exit('\nUSAGE ERROR: Please enter either yes or no')
        ###Si/Ge heterostructures
        if split_region == 'yes':
            left_ids = 0
            lint_ids= 0
            for i in range(len(intdata)):
                if intdata[i,1] == 2:
                    lint_ids = i
            lint_data = intdata[:lint_ids+1,:]
            intdata = np.delete(intdata, range(lint_ids+1),axis=0)
            rint_data = intdata
            del intdata
            for i in range(len(region)):
                if region[i][1] == 2:
                    left_ids = i
            left_region = region[:left_ids+1][:]
            si_atoms = len(left_region)
            region = np.delete(region, range(left_ids+1), axis=0)
            right_region = region
            ge_atoms = len(right_region)
            del region
            print ('\n\tThere are ' + str(len(left_region)) + ' atoms in the Si section and ' +
                   str(len(right_region)) + ' atoms in the Ge region.\n')
            ###Si side
            try:
                ###find out how many point defects for heterogeneous systems
                l_type = raw_input('\nWhat is the species of the Si region subsitutional? (USAGE: none or C, Si, or Ge)\n>')
                if l_type != 'C' and l_type != 'Si' and l_type != 'Ge' and l_type != 'none':
                    sys.exit('USAGE ERROR: Substitutional type must be "none" or C, Si, or Ge')
                lnum_subs = int(raw_input('\nHow many substitutional defects in Si region?\n>'))
                lnum_self = int(raw_input('\nHow many SELF-interstitials in Si region?\n>'))
                lnum_ints = int(raw_input('\nHow many IMPURITY-interstitials in Si region?\n'
                                          'NOTE: impurity type is same as subsitutional.\n>'))
                lnum_vacs = int(raw_input('\nHow many vacancies in Si region?\n>'))
                ###find out if pairs are included
                lnum_selffrenkel = int(raw_input('\nHow many SELF-intersitial-vacancy pairs in Si region?\n>'))
                lnum_frenkel = int(raw_input('\nHow many IMPURTIY-interstital-vacancy pairs in Si region?\n'
                                             'NOTE: impurity type is same as subsitutional.\n>'))
                lnum_v2 = int(raw_input('\nHow many vacancy pairs in Si region?\n>'))
            except ValueError:
                sys.exit('\nEnter an integer between 0 and ' + str(len(left_region)) + '\n')
            ###Ge side
            try:
                ###find out how many point defects for heterogeneous systems
                r_type = raw_input('\n---------------------------------------------------------------------------------\n'
                                   '\nWhat is the species of the Ge region subsitutional? (USAGE: none or C, Si, or Ge)\n>')
                if r_type != 'C' and r_type != 'Si' and r_type != 'Ge' and r_type != 'none':
                    sys.exit('USAGE ERROR: Substitutional type must be "none" or C, Si, or Ge')
                rnum_subs = int(raw_input('\nHow many substitutional defects in Ge region?\n>'))
                rnum_self = int(raw_input('\nHow many SELF-interstitials in Ge region?\n>'))
                rnum_ints = int(raw_input('\nHow many IMPURITY-interstitials in Ge region?\n'
                                          'NOTE: impurity type is same as subsitutional.\n>'))
                rnum_vacs = int(raw_input('\nHow many vacancies in Ge region?\n>'))
                ###find out if pairs are included
                rnum_selffrenkel = int(raw_input('\nHow many SELF-intersitial-vacancy pairs in Ge region?\n>'))
                rnum_frenkel = int(raw_input('\nHow many IMPURTIY-interstital-vacancy pairs in Ge region?\n'
                                             'NOTE: impurity type is same as subsitutional.\n>'))
                rnum_v2 = int(raw_input('\nHow many vacancy pairs in Ge region?\n>'))
            except ValueError:
                sys.exit('\nEnter an integer between 0 and ' + str(len(right_region)) + '\n')
            
    ###homogeneous systems
    if split_region == 'no':
        print ('\nThere are ' + str(len(region)) + ' atoms in the defected region.\n')
        try:
            ###find out how many point defects for homogeneous systems
            all_type = raw_input('\nWhat is the species of the subsitutional? (USAGE: none or C, Si, or Ge)\n>')
            if all_type != 'C' and all_type != 'Si' and all_type != 'Ge' and all_type != 'none':
                    sys.exit('USAGE ERROR: Substitutional type must be "none" or C, Si, or Ge')
            num_subs = int(raw_input('\nHow many substitutional defects?\n>'))
            num_self = int(raw_input('\nHow many SELF-interstitials?\n>'))
            num_ints = int(raw_input('\nHow many IMPURITY-interstitials?\n'
                                     'NOTE: impurity type is same as subsitutional.\n>'))
            num_vacs = int(raw_input('\nHow many vacancies?\n>'))
            ###find out if pairs are included
            num_selffrenkel = int(raw_input('\nHow many SELF-intersitial-vacancy pairs?\n>'))
            num_frenkel = int(raw_input('\nHow many IMPURTIY-interstital-vacancy pairs?\n'
                                        'NOTE: impurity type is same as subsitutional.\n>'))
            num_v2 = int(raw_input('\nHow many vacancy pairs?\n>'))
        except ValueError:
            sys.exit('\nEnter an integer between 0 and ' + len(region) + '\n')
####################################################





########    PLACE DEFECTS    #######################
if str(defect) != 'no':
    example = cp.deepcopy(data)
    ###interfacial mixing
    if split_region == 'yes': #if specfied to have different defects on each side
        
        ###left side
        if lnum_subs != 0: #remove subs from region and change type
            randids = np.array(random.sample(range(len(left_region)), lnum_subs))
            lsubs = left_region[randids,:]
            left_region = np.delete(left_region, randids, axis=0)
            lsubs[:,1] = 8
            data = np.append(data,lsubs,axis=0)

        if lnum_self != 0: #remove self interstitial sites from region
            randids = np.array(random.sample(range(len(lint_data)), lnum_self))
            lself = lint_data[randids,:]
            lint_data = np.delete(lint_data,randids,axis=0)
            lself[:,1] = left_region[0,1]
            data = np.append(data,lself,axis=0)
        
        if lnum_ints != 0: #remove impurity interstitial sites from region
            randids = np.array(random.sample(range(len(lint_data)), lnum_ints))
            lints = lint_data[randids,:]
            lint_data = np.delete(lint_data,randids,axis=0)
            lints[:,1] = 8
            data = np.append(data,lints,axis=0)
        
        if lnum_vacs != 0: #remove single vacancy sites from region
            randids = np.array(random.sample(range(len(left_region)), lnum_vacs))
            lvacs = left_region[randids,:]
            left_region = np.delete(left_region, randids, axis=0)
            lvacs[:,1] = 10
            example = np.append(example,lvacs,axis=0)
            
        if lnum_selffrenkel != 0: #remove self-inter sites to form pairs
            randids = np.array(random.sample(range(len(lint_data)), lnum_selffrenkel))
            lselffrenkel = lint_data[randids,:]
            lint_data = np.delete(lint_data,randids,axis=0)
            lselffrenkel[:,1] = left_region[0,1]
            lselffrenkelpair = np.zeros((len(lselffrenkel),5))
            for i in range(len(lselffrenkel)):
                nn_dist = np.zeros((len(left_region),2))
                for j in range(len(left_region)):
                    nn_dist[j,0] = j
                    nn_dist[j,1] = math.sqrt((lselffrenkel[i,2]-left_region[j,2])**2+(lselffrenkel[i,3]-left_region[j,3])**2
                                   +(lselffrenkel[i,4]-left_region[j,4])**2)
                nn_dist = nn_dist[nn_dist[:,1].argsort()]
                num = (nn_dist[:,1] <= nn_dist[0,1]+0.1).sum()
                rand = random.sample(range(num),1)
                lselffrenkelpair[i,:] = left_region[int(nn_dist[rand,0]),:] #vacancy pair
                left_region = np.delete(left_region, int(nn_dist[rand,0]), axis=0) #delete from population
            lselffrenkelpair[:,1] = 10
            data = np.append(data,lselffrenkel,axis=0)
            example = np.append(example,lselffrenkelpair,axis=0)
        
        if lnum_frenkel != 0: #remove inter sites to form pairs
            randids = np.array(random.sample(range(len(lint_data)), lnum_frenkel))
            lfrenkel = lint_data[randids,:]
            lint_data = np.delete(lint_data,randids,axis=0)
            lfrenkel[:,1] = 8
            lfrenkelpair = np.zeros((len(lfrenkel),5))
            for i in range(len(lfrenkel)):
                nn_dist = np.zeros((len(left_region),2))
                for j in range(len(left_region)):
                    nn_dist[j,0] = j
                    nn_dist[j,1] = math.sqrt((lfrenkel[i,2]-left_region[j,2])**2+(lfrenkel[i,3]-left_region[j,3])**2
                                   +(lfrenkel[i,4]-left_region[j,4])**2)
                nn_dist = nn_dist[nn_dist[:,1].argsort()]
                num = (nn_dist[:,1] <= nn_dist[0,1]+0.1).sum()
                rand = random.sample(range(num),1)
                lfrenkelpair[i,:] = left_region[int(nn_dist[rand,0]),:] #vacancy pair
                left_region = np.delete(left_region, int(nn_dist[rand,0]), axis=0) #delete from population
            lfrenkelpair[:,1] = 10
            data = np.append(data,lfrenkel,axis=0)
            example = np.append(example,lfrenkelpair,axis=0)
        
        if lnum_v2 != 0: #remove vacancy sites to form pairs
            randids = np.array(random.sample(range(len(left_region)), lnum_v2))
            lv2 = left_region[randids,:]
            left_region = np.delete(left_region, randids, axis=0)
            lv2pair = np.zeros((len(lv2),5))
            for i in range(len(lv2)):
                nn_dist = np.zeros((len(left_region),2))
                for j in range(len(left_region)):
                    nn_dist[j,0] = j
                    nn_dist[j,1] = math.sqrt((lv2[i,2]-left_region[j,2])**2+(lv2[i,3]-left_region[j,3])**2
                                   +(lv2[i,4]-left_region[j,4])**2)
                nn_dist = nn_dist[nn_dist[:,1].argsort()]
                num = (nn_dist[:,1] <= nn_dist[0,1]+0.1).sum()
                rand = random.sample(range(num),1)
                lv2pair[i,:] = left_region[int(nn_dist[rand,0]),:] #second vacancies
                left_region = np.delete(left_region, int(nn_dist[rand,0]), axis=0) #delete from population
            lv2[:,1] = 10
            lv2pair[:,1] = 10
            example = np.append(example,lv2,axis=0)
            example = np.append(example,lv2pair,axis=0)

        
        ###right side
        if rnum_subs != 0:
            randids = np.array(random.sample(range(len(right_region)), rnum_subs))
            rsubs = right_region[randids,:]
            right_region = np.delete(right_region, randids, axis=0)
            rsubs[:,1] = 7
            data = np.append(data,rsubs,axis=0)
        
        if rnum_self != 0:
            randids = np.array(random.sample(range(len(rint_data)), rnum_self))
            rself = rint_data[randids,:]
            rint_data = np.delete(rint_data,randids,axis=0)
            rself[:,1] = right_region[0,1]
            data = np.append(data,rself,axis=0)
        
        if rnum_ints != 0:
            randids = np.array(random.sample(range(len(rint_data)), rnum_ints))
            rints = rint_data[randids,:]
            rint_data = np.delete(rint_data,randids,axis=0)
            rints[:,1] = 7
            data = np.append(data,rints,axis=0)
        
        if rnum_vacs != 0:
            randids = np.array(random.sample(range(len(right_region)), rnum_vacs))
            rvacs = right_region[randids,:]
            right_region = np.delete(right_region, randids, axis=0)
            rvacs[:,1] = 10
            example = np.append(example,rvacs,axis=0)
        
        if rnum_selffrenkel != 0:
            randids = np.array(random.sample(range(len(rint_data)), rnum_selffrenkel))
            rselffrenkel = rint_data[randids,:]
            rint_data = np.delete(rint_data,randids,axis=0)
            rselffrenkel[:,1] = right_region[0,1]
            rselffrenkelpair = np.zeros((len(rselffrenkel),5))
            for i in range(len(rselffrenkel)):
                nn_dist = np.zeros((len(right_region),2))
                for j in range(len(right_region)):
                    nn_dist[j,0] = j
                    nn_dist[j,1] = math.sqrt((rselffrenkel[i,2]-right_region[j,2])**2+(rselffrenkel[i,3]-right_region[j,3])**2
                                   +(rselffrenkel[i,4]-right_region[j,4])**2)
                nn_dist = nn_dist[nn_dist[:,1].argsort()]
                num = (nn_dist[:,1] <= nn_dist[0,1]+0.1).sum()
                rand = random.sample(range(num),1)
                rselffrenkelpair[i,:] = right_region[int(nn_dist[rand,0]),:] #vacancy pair
                right_region = np.delete(right_region, int(nn_dist[rand,0]), axis=0) #delete from population
            rselffrenkelpair[:,1] = 10
            data = np.append(data,rselffrenkel,axis=0)
            example = np.append(example,rselffrenkelpair,axis=0)
        
        if rnum_frenkel != 0:
            randids = np.array(random.sample(range(len(rint_data)), rnum_frenkel))
            rfrenkel = rint_data[randids,:]
            rint_data = np.delete(rint_data,randids,axis=0)
            rfrenkel[:,1] = 7
            rfrenkelpair = np.zeros((len(rfrenkel),5))
            for i in range(len(rfrenkel)):
                nn_dist = np.zeros((len(right_region),2))
                for j in range(len(right_region)):
                    nn_dist[j,0] = j
                    nn_dist[j,1] = math.sqrt((rfrenkel[i,2]-right_region[j,2])**2+(rfrenkel[i,3]-right_region[j,3])**2
                                   +(rfrenkel[i,4]-right_region[j,4])**2)
                nn_dist = nn_dist[nn_dist[:,1].argsort()]
                num = (nn_dist[:,1] <= nn_dist[0,1]+0.1).sum()
                rand = random.sample(range(num),1)
                rfrenkelpair[i,:] = right_region[int(nn_dist[rand,0]),:] #vacancy pair
                right_region = np.delete(right_region, int(nn_dist[rand,0]), axis=0) #delete from population
            rfrenkelpair[:,1] = 10
            data = np.append(data,rfrenkel,axis=0)
            example = np.append(example,rfrenkelpair,axis=0)

        if rnum_v2 != 0:
            randids = np.array(random.sample(range(len(right_region)), rnum_v2))
            rv2 = right_region[randids,:]
            right_region = np.delete(right_region, randids, axis=0)
            rv2pair = np.zeros((len(rv2),5))
            for i in range(len(rv2)):
                nn_dist = np.zeros((len(right_region),2))
                for j in range(len(right_region)):
                    nn_dist[j,0] = j
                    nn_dist[j,1] = math.sqrt((rv2[i,2]-right_region[j,2])**2+(rv2[i,3]-right_region[j,3])**2
                                   +(rv2[i,4]-right_region[j,4])**2)
                nn_dist = nn_dist[nn_dist[:,1].argsort()]
                num = (nn_dist[:,1] <= nn_dist[0,1]+0.1).sum()
                rand = random.sample(range(num),1)
                rv2pair[i,:] = right_region[int(nn_dist[rand,0]),:] #second vacancies
                right_region = np.delete(right_region, int(nn_dist[rand,0]), axis=0) #delete from population
            rv2[:,1] = 10
            rv2pair[:,1] = 10
            example = np.append(example,rv2,axis=0)
            example = np.append(example,rv2pair,axis=0)
            
        data = np.append(np.append(data,left_region,axis=0),right_region,axis=0)
      
        
    ###uniform defects
    else:
        if num_subs != 0:
            randids = np.array(random.sample(range(len(region)), num_subs))
            allsubs = region[randids,:]
            region = np.delete(region, randids, axis=0)
            allsubs[:,1] = 9
            data = np.append(data,allsubs,axis=0)
                
        if num_self != 0:
            randids = np.array(random.sample(range(len(intdata)), num_self))
            allself = intdata[randids,:]
            intdata = np.delete(intdata,randids,axis=0)
            allself[:,1] = region[0,1]
            data = np.append(data,allself,axis=0)
        
        if num_ints != 0:
            randids = np.array(random.sample(range(len(intdata)), num_ints))
            allints = intdata[randids,:]
            intdata = np.delete(intdata,randids,axis=0)
            allints[:,1] = 9
            data = np.append(data,allints,axis=0)
        
        if num_vacs != 0:
            randids = np.array(random.sample(range(len(region)), num_vacs))
            allvacs = region[randids,:]
            region = np.delete(region, randids, axis=0)
            allvacs[:,1] = 10
            example = np.append(example,allvacs,axis=0)
        
        if num_selffrenkel != 0:
            randids = np.array(random.sample(range(len(intdata)), num_selffrenkel))
            allselffrenkel = intdata[randids,:]
            intdata = np.delete(intdata,randids,axis=0)
            allselffrenkel[:,1] = region[0,1]
            allselffrenkelpair = np.zeros((len(allselffrenkel),5))
            for i in range(len(allselffrenkel)):
                nn_dist = np.zeros((len(region),2))
                for j in range(len(region)):
                    nn_dist[j,0] = j
                    nn_dist[j,1] = math.sqrt((allselffrenkel[i,2]-region[j,2])**2+(allselffrenkel[i,3]-region[j,3])**2
                                   +(allselffrenkel[i,4]-region[j,4])**2)
                nn_dist = nn_dist[nn_dist[:,1].argsort()]
                num = (nn_dist[:,1] <= nn_dist[0,1]+0.1).sum()
                rand = random.sample(range(num),1)
                allselffrenkelpair[i,:] = region[int(nn_dist[rand,0]),:] #vacancy pair
                region = np.delete(region, int(nn_dist[rand,0]), axis=0) #delete from population
            allselffrenkelpair[:,1] = 10
            data = np.append(data,allselffrenkel,axis=0)
            example = np.append(example,allselffrenkelpair,axis=0)
        
        if num_frenkel != 0:
            randids = np.array(random.sample(range(len(intdata)), num_frenkel))
            allfrenkel = intdata[randids,:]
            intdata = np.delete(intdata,randids,axis=0)
            allfrenkel[:,1] = 9
            allfrenkelpair = np.zeros((len(allfrenkel),5))
            for i in range(len(allfrenkel)):
                nn_dist = np.zeros((len(region),2))
                for j in range(len(region)):
                    nn_dist[j,0] = j
                    nn_dist[j,1] = math.sqrt((allfrenkel[i,2]-region[j,2])**2+(allfrenkel[i,3]-region[j,3])**2
                                   +(allfrenkel[i,4]-region[j,4])**2)
                nn_dist = nn_dist[nn_dist[:,1].argsort()]
                num = (nn_dist[:,1] <= nn_dist[0,1]+0.1).sum()
                rand = random.sample(range(num),1)
                allfrenkelpair[i,:] = region[int(nn_dist[rand,0]),:] #vacancy pair
                region = np.delete(region, int(nn_dist[rand,0]), axis=0) #delete from population
            allfrenkelpair[:,1] = 10
            data = np.append(data,allfrenkel,axis=0)
            example = np.append(example,allfrenkelpair,axis=0)
                
        if num_v2 != 0:
            randids = np.array(random.sample(range(len(region)), num_v2))
            allv2 = region[randids,:]
            region = np.delete(region, randids, axis=0)
            allv2pair = np.zeros((len(allv2),5))
            for i in range(len(allv2)):
                nn_dist = np.zeros((len(region),2))
                for j in range(len(region)):
                    nn_dist[j,0] = j
                    nn_dist[j,1] = math.sqrt((allv2[i,2]-region[j,2])**2+(allv2[i,3]-region[j,3])**2
                                   +(allv2[i,4]-region[j,4])**2)
                nn_dist = nn_dist[nn_dist[:,1].argsort()]
                num = (nn_dist[:,1] <= nn_dist[0,1]+0.1).sum()
                rand = random.sample(range(num),1)
                allv2pair[i,:] = region[int(nn_dist[rand,0]),:] #second vacancies
                region = np.delete(region, int(nn_dist[rand,0]), axis=0) #delete from population
            allv2[:,1] = 10
            allv2pair[:,1] = 10
            example = np.append(example,allv2,axis=0)
            example = np.append(example,allv2pair,axis=0)     
        
        data = np.append(data,region,axis=0)

    example = np.unique(np.append(example,data,axis=0),axis=0)
####################################################




            
########    SCALE LATTICE AND PLOT    ##############
###apply lattice constant        
data = data[np.lexsort((data[:,4],data[:,3],data[:,2]))]

for i in range(len(data)):
    data[i][2] = data[i][2]*b/a
    data[i][3] = data[i][3]*b/a
    data[i][4] = data[i][4]*b/a
     
if defect == 'yes':     
    for i in range(len(example)):
        example[i][2] = example[i][2]*b/a
        example[i][3] = example[i][3]*b/a
        example[i][4] = example[i][4]*b/a
          
###find max and min coords          
maxcoords = data.max(axis=0)
mincoords = data.min(axis=0)
xmax = maxcoords[2]
ymax = maxcoords[3]
zmax = maxcoords[4]
xmin = mincoords[2]
ymin = mincoords[3]
zmin = mincoords[4]
del maxcoords, mincoords
#####################################################





###########  WRITE CONFIG INFO TO TXT FILE  #########
with open('config_info.txt', 'w') as f:
    if defect == 'yes':
        f.write('Configuration information for FCC-diamond structure of type "' + species + '"\n\n')
        f.write('Original number of atoms was ' + str(size) + '\n')
        f.write('Total number of remaining atoms is ' + str(len(data)) + '\n\n')
        f.write('Defect concentrations are given as percent of the SPECIFIED REGION, i.e atoms excluded from the defected region are not counted.\n')
        f.write('TOTAL defect concentration is listed at the END.\n\n')
        if split_region == 'yes':
            atoms = si_atoms+ge_atoms
            total = (lnum_subs+lnum_self+lnum_ints+lnum_vacs+lnum_selffrenkel*2+lnum_frenkel*2+lnum_v2*2+
                     rnum_subs+rnum_self+rnum_ints+rnum_vacs+rnum_selffrenkel*2+rnum_frenkel*2+rnum_v2*2)
            f.write('There were ' + str(si_atoms) + ' Si atoms and ' + str(ge_atoms) + ' Ge atoms selected'
                    ' in which to place defects\n\n')
            f.write('Si region defects:\n')
            f.write('\tIMPURITY TYPE =\t\t\t' + l_type + '\n')
            f.write('\tSUBSTITUTIONALS =\t\t' + str(lnum_subs) + ' =\t\t' + str(100.0*lnum_subs/atoms) + '%\n')
            f.write('\tSELF-INTERSTITIALS =\t\t' + str(lnum_self) + ' =\t\t' + str(100.0*lnum_self/atoms) + '%\n')
            f.write('\tIMPURITY-INTERSTITIALS =\t' + str(lnum_ints) + ' =\t\t' + str(100.0*lnum_ints/atoms) + '%\n')
            f.write('\tVACANCIES =\t\t\t' + str(lnum_vacs) + ' =\t\t' + str(100.0*lnum_vacs/atoms) + '%\n')
            f.write('\tSELF-INTERSTITIAL-PAIRS =\t' + str(lnum_selffrenkel) + ' =\t\t' + str(100.0*2*lnum_selffrenkel/atoms) + '% (x2 for pairs)\n')
            f.write('\tIMPURITY-INTERSTITIAL-PAIRS =\t' + str(lnum_frenkel) + ' =\t\t' + str(100.0*2*lnum_frenkel/atoms) + ' % (x2 for pairs)\n')
            f.write('\tVACANCY-PAIRS =\t\t\t' + str(lnum_v2) + ' =\t\t' + str(100.0*2*lnum_v2/atoms) + '% (x2 for pairs)\n\n')
    
            f.write('Ge region defects:\n')
            f.write('\tIMPURITY TYPE =\t\t\t' + r_type + '\n')
            f.write('\tSUBSTITUTIONALS =\t\t' + str(rnum_subs) + ' =\t\t' + str(100.0*rnum_subs/atoms) + '%\n')
            f.write('\tSELF-INTERSTITIALS =\t\t' + str(rnum_self) + ' =\t\t' + str(100.0*rnum_self/atoms) + '%\n')
            f.write('\tIMPURITY-INTERSTITIALS =\t' + str(rnum_ints) + ' =\t\t' + str(100.0*rnum_ints/atoms) + '%\n')
            f.write('\tVACANCIES =\t\t\t' + str(rnum_vacs) + ' =\t\t' + str(100.0*rnum_vacs/atoms) + '%\n')
            f.write('\tSELF-INTERSTITIAL-PAIRS =\t' + str(rnum_selffrenkel) + ' =\t\t' + str(100.0*2*rnum_selffrenkel/atoms) + '% (x2 for pairs)\n')
            f.write('\tIMPURITY-INTERSTITIAL-PAIRS =\t' + str(rnum_frenkel) + ' =\t\t' + str(100.0*2*rnum_frenkel/atoms) + ' % (x2 for pairs)\n')
            f.write('\tVACANCY-PAIRS =\t\t\t' + str(rnum_v2) + ' =\t\t' + str(100.0*2*rnum_v2/atoms) + '% (x2 for pairs)\n\n')
            f.write('\tPARTIAL-DEFECTS =\t\t' + str(total) + '/' + str(atoms) + ' =\t' + str(100.0*total/atoms) + '%\n')
            f.write('\tTOTAL-DEFECTS =\t\t\t' + str(total) + '/' + str(size) + ' =\t' + str(100.0*total/size) + '%')
            
        else:
            total = (num_subs+num_self+num_ints+num_vacs+num_selffrenkel*2+num_frenkel*2+num_v2*2)
            f.write('There were ' + str(atoms) + ' ' + species + ' atoms selected'
                    ' in which to place defects\n\n')
            f.write('Defects:\n')
            f.write('\tIMPURITY TYPE =\t\t\t' + all_type + '\n')
            f.write('\tSUBSTITUTIONALS =\t\t' + str(num_subs) + ' =\t\t' + str(100.0*num_subs/atoms) + '%\n')
            f.write('\tSELF-INTERSTITIALS =\t\t' + str(num_self) + ' =\t\t' + str(100.0*num_self/atoms) + '%\n')
            f.write('\tIMPURITY-INTERSTITIALS =\t' + str(num_ints) + ' =\t\t' + str(100.0*num_ints/atoms) + '%\n')
            f.write('\tVACANCIES =\t\t\t' + str(num_vacs) + ' =\t\t' + str(100.0*num_vacs/atoms) + '%\n')
            f.write('\tSELF-INTERSTITIAL-PAIRS =\t' + str(num_selffrenkel) + ' =\t\t' + str(100.0*2*num_selffrenkel/atoms) + '% (x2 for pairs)\n')
            f.write('\tIMPURITY-INTERSTITIAL-PAIRS =\t' + str(num_frenkel) + ' =\t\t' + str(100.0*2*num_frenkel/atoms) + ' % (x2 for pairs)\n')
            f.write('\tVACANCY-PAIRS =\t\t\t' + str(num_v2) + ' =\t\t' + str(100.0*2*num_v2/atoms) + '% (x2 for pairs)\n\n')
            f.write('\tPARTIAL-DEFECTS =\t\t' + str(total) + '/' + str(atoms) + ' =\t' + str(100.0*total/atoms) + '%\n')
            f.write('\tTOTAL-DEFECTS =\t\t\t' + str(total) + '/' + str(size) + ' =\t' + str(100.0*total/size) + '%')
            
    else:
        f.write('NO DEFECTS')
#####################################################





###########  WRITE TO .xyz FILE    ##################
if defect == 'yes':
    if (species == 'Si/Ge'):
        exfilename = ('EXAMPLE.Si.' + str(n_si) + '.Ge.' + str(n_ge) + 'x' + str(ny) + 'x' +str(nz) + '.xyz')
    else:
        exfilename = ('EXAMPLE.' + species + '.' + str(nx) + 'x' + str(ny) + 'x' +str(nz) + '.xyz')
    with open(exfilename, 'w') as f:
        f.write(str(len(example))+'\n')
        f.write(str(xmax)+'\t0'+'\t0'+'\t0\t'+str(ymax)+'\t0'+'\t0'+'\t0\t'+str(zmax)+'\n')
        for i in range(len(example)-1):
            if example[i][1] == 1:
                f.write('C\t'+str(example[i][2])+'\t'+str(example[i][3])+'\t'+str(example[i][4])+'\n')
            if example[i][1] == 2:
                f.write('Si\t'+str(example[i][2])+'\t'+str(example[i][3])+'\t'+str(example[i][4])+'\n')
            if example[i][1] == 3:
                f.write('Ge\t'+str(example[i][2])+'\t'+str(example[i][3])+'\t'+str(example[i][4])+'\n')
            if example[i][1] == 9:
                f.write(all_type+'\t'+str(example[i][2])+'\t'+str(example[i][3])+'\t'+str(example[i][4])+'\n')
            if example[i][1] == 8:
                f.write(l_type+'\t'+str(example[i][2])+'\t'+str(example[i][3])+'\t'+str(example[i][4])+'\n')
            if example[i][1] == 7:
                f.write(r_type+'\t'+str(example[i][2])+'\t'+str(example[i][3])+'\t'+str(example[i][4])+'\n')
            if example[i][1] == 10:
                f.write('Li\t' + str(example[i][2]) + '\t' +str(example[i][3]) + '\t' + str(example[i][4]) + '\n')
        if example[-1][1] == 1:
            f.write('C\t'+str(example[-1][2])+'\t'+str(example[-1][3])+'\t'+str(example[-1][4]))
        if example[-1][1] == 2:
            f.write('Si\t'+str(example[-1][2])+'\t'+str(example[-1][3])+'\t'+str(example[-1][4]))
        if example[-1][1] == 3:
            f.write('Ge\t'+str(example[-1][2])+'\t'+str(example[-1][3])+'\t'+str(example[-1][4]))
        if example[-1][1] == 9:
            f.write(all_type+'\t'+str(example[-1][2])+'\t'+str(example[-1][3])+'\t'+str(example[-1][4]))
        if example[-1][1] == 8:
            f.write(l_type+'\t'+str(example[-1][2])+'\t'+str(example[-1][3])+'\t'+str(example[-1][4]))
        if example[-1][1] == 7:
            f.write(r_type+'\t'+str(example[-1][2])+'\t'+str(example[-1][3])+'\t'+str(example[-1][4]))
        if example[-1][1] == 10:
            f.write('Li\t' + str(example[-1][2]) + '\t' +str(example[-1][3]) + '\t' + str(example[-1][4]))


if (species == 'Si/Ge'):
    filename = ('Si.' + str(n_si) + '.Ge.' + str(n_ge) + 'x' + str(ny) + 'x' +str(nz) + '.xyz')
else:
    filename = (species + '.' + str(nx) + 'x' + str(ny) + 'x' +str(nz) + '.xyz')
with open(filename, 'w') as f:
    f.write(str(len(data))+'\n')
    f.write(str(xmax)+'\t0'+'\t0'+'\t0\t'+str(ymax)+'\t0'+'\t0'+'\t0\t'+str(zmax)+'\n')
    for i in range(len(data)-1):
        if data[i][1] == 1:
            f.write('C\t'+str(data[i][2])+'\t'+str(data[i][3])+'\t'+str(data[i][4])+'\n')
        if data[i][1] == 2:
            f.write('Si\t'+str(data[i][2])+'\t'+str(data[i][3])+'\t'+str(data[i][4])+'\n')
        if data[i][1] == 3:
            f.write('Ge\t'+str(data[i][2])+'\t'+str(data[i][3])+'\t'+str(data[i][4])+'\n')
        if data[i][1] == 9:
            f.write(all_type+'\t'+str(data[i][2])+'\t'+str(data[i][3])+'\t'+str(data[i][4])+'\n')
        if data[i][1] == 8:
            f.write(l_type+'\t'+str(data[i][2])+'\t'+str(data[i][3])+'\t'+str(data[i][4])+'\n')
        if data[i][1] == 7:
            f.write(r_type+'\t'+str(data[i][2])+'\t'+str(data[i][3])+'\t'+str(data[i][4])+'\n')
    if data[-1][1] == 1:
        f.write('C\t'+str(data[-1][2])+'\t'+str(data[-1][3])+'\t'+str(data[-1][4]))
    if data[-1][1] == 2:
        f.write('Si\t'+str(data[-1][2])+'\t'+str(data[-1][3])+'\t'+str(data[-1][4]))
    if data[-1][1] == 3:
        f.write('Ge\t'+str(data[-1][2])+'\t'+str(data[-1][3])+'\t'+str(data[-1][4]))
    if data[-1][1] == 9:
        f.write(all_type+'\t'+str(data[-1][2])+'\t'+str(data[-1][3])+'\t'+str(data[-1][4]))
    if data[-1][1] == 8:
        f.write(l_type+'\t'+str(data[-1][2])+'\t'+str(data[-1][3])+'\t'+str(data[-1][4]))
    if data[-1][1] == 7:
        f.write(r_type+'\t'+str(data[-1][2])+'\t'+str(data[-1][3])+'\t'+str(data[-1][4]))
####################################################




###########  WRITE LAMMPS DATA  ###################
buff = b/8.0
if species == 'Si/Ge':
    datafile = ('data.Si.' + str(n_si) + '.Ge.' + str(n_ge) + 'x' + str(ny) + 'x' +str(nz))
else:
    datafile = ('data.' + species + '.' + str(nx) + 'x' + str(ny) + 'x' +str(nz))

for j in range(len(data)):
    if data[j,1] == 9:
        if all_type == 'C':
            data[j,1] = 1
        elif all_type == 'Si':
            data[j,1] = 2
        elif all_type == 'Ge':
            data[j,1] = 3
    elif data[j,1] == 8:
        if l_type == 'C':
            data[j,1] = 1
        elif l_type == 'Si':
            data[j,1] = 2
        elif l_type == 'Ge':
            data[j,1] = 3
    elif data[j,1] == 7:
        if r_type == 'C':
            data[j,1] = 1
        elif r_type == 'Si':
            data[j,1] = 2
        elif r_type == 'Ge':
            data[j,1] = 3

types = np.zeros((len(data),1))
types = data[:,1]
types = np.unique(types,axis=0)
    
masses = np.zeros((len(types),1))
for i in range(len(masses)):
    if types[i] == 1:
        masses[i] = 12.0107 #mass of C
    elif types[i] == 2:
        masses[i] = 28.0855 #mass of si
    elif types[i] == 3:
        masses[i] = 72.6400 #mass of ge
 
for i in range(len(types)):
    for j in range(len(data)):
        if data[j,1] == types[i]:
            data[j,1] = -(i+1)
            
for i in range(len(data)):
    data[i,1] = abs(data[i,1])
    
    
with open(datafile, 'w') as f:
    
    f.write(str('LAMMPS DATA FILE\n'))
    
    f.write('\n' + str(len(data)) + ' atoms\n')
    f.write('\n' + str(len(masses)) + ' atom types\n')
    f.write('\n' + str(xmin-buff)+' '+str(xmax+buff)+' xlo'+' xhi\n')
    f.write(str(ymin-buff)+' '+str(ymax+buff)+' ylo'+' yhi\n')
    f.write(str(zmin-buff)+' '+str(zmax+buff)+' zlo'+' zhi\n')
    f.write('\nMasses\n')
    for i in range(len(masses)):
        f.write('\n' + str(i+1) + ' ' + str(float(masses[i])))
    f.write('\n\nAtoms\n\n')
    for i in range(len(data)-1):
        f.write(str(int(i+1)) + ' ' + str(int(data[i,1])) + ' ' + str(data[i,2]) + ' ' +
                str(data[i,3]) + ' ' + str(data[i,4]) + '\n')
    f.write(str(len(data)) +  ' ' + str(int(data[-1,1])) + ' ' + str(data[-1,2]) + ' ' + 
            str(data[-1,3]) + ' ' + str(data[-1,4]))
##################################################

if defect != 'yes':     
    print('\n\t-------------------------------------------------------\n\t'
          'ALL DONE! Find "' + filename + '", "config_info.txt",\n\t'
          'and "' + datafile + '" in your current directory!\n\n\t"config_info.txt" contains information about'
          ' the numbers and types of defects in your system.\n\tKeep it with the config files!\n\n\n')        
else:
    print('\n\t-------------------------------------------------------\n\t'
          'ALL DONE! Find "' + filename + '", "' + exfilename + '", "config_info.txt",\n\t'
          'and "' + datafile + '" in your current directory!\n\n\t"config_info.txt" contains information about'
          ' the numbers and types of defects in your system.\n\tKeep it with the config files!\n\n\n')        

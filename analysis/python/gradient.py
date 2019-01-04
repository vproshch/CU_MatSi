#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  6 15:24:38 2018

@author: ty
"""

import numpy as np
import sys
import matplotlib.pyplot as plt
#plt.switch_backend('TKagg')

##### GET TEMP DATA #####
suffix = raw_input("Enter the suffix of the temperature file:\n->")
filename = 'tmp.profile.'+suffix
num_lines = sum(1 for line in open(filename))-3 #how many lines

with open(filename, 'r') as fid:
        
    for i in range(3): #ignore comments
        fid.readline()
        
    info = fid.readline().strip().split() #get number of chunks
    num_blocks = num_lines/(int(info[1])+1) #get number of temp. blocks
        
        
    data = np.zeros((int(info[1]),3)) #initialze array for temp. data
    
    #only keep the last 10% of the data... converged
    num_keep = np.round(int(num_blocks*0.9))
    for i in range(num_blocks-num_keep-1):
        for j in range(int(info[1])+1):
            etc = fid.readline()
    for i in range(int(info[1])):
        etc = fid.readline()
        
    for i in range(num_keep):
        etc = fid.readline()
        for j in range(int(info[1])):
            tmp = fid.readline().strip().split()
            if i == 0:
                data[j,0] = int(tmp[0])
                data[j,1] = float(tmp[1])
            data[j,2] = data[j,2]+float(tmp[3])
    data[:,2] = data[:,2]/num_keep
##############
    
##### PLOT TMP PROFILE #####
print('\nEnter chunk indicies from plot for fitting.'
      '\nUsage: (Enter 5 points) 1: The left side of the left chunk, 2: the \n'
      'right side of the left chunk, 3: the center of the interface if present,\n'
      '4: the left side of the right chunk, and 5: the right side of the right chunk.\n')
    
plt.plot(data[:,0],data[:,2],'kx--')
plt.xlabel('Chunk index')
plt.ylabel('Temperature, K')
plt.grid(b=None, which='major', axis='both')
plt.show()
############################

##### GET INDICIES AND FIT #####
ind = np.zeros((5,1))
ind[0] = int(raw_input("\nLeft side of left chunk:\n->"))
if ind[0] < 0:
    sys.exit('Index must be positive')
ind[1] = int(raw_input("\nRight side of left chunk:\n->"))
if ind[1] <= ind[0]:
    sys.exit('Index must greater than previous one')
ind[2] = int(raw_input("\nCenter of interface (if present):\n->"))
if ind[2] <= ind[1]:
    sys.exit('Index must greater than previous one')
ind[3] = int(raw_input("\nLeft side of right chunk:\n->"))
if ind[3] <= ind[2]:
    sys.exit('Index must greater than previous one')
ind[4] = int(raw_input("\nRight side of right chunk:\n->"))
if ind[4] <= ind[3]:
    sys.exit('Index must greater than previous one')
    
left = data[(int(ind[0])-1):int(ind[1]),:]
right = data[(int(ind[3])-1):int(ind[4]),:]

fit_left = np.polyfit(left[:,1],left[:,2],1)
fit_right = np.polyfit(right[:,1],right[:,2],1)

stretch_left = np.polyval(fit_left,data[int(ind[0]-1):int(ind[2]),1])
stretch_right = np.polyval(fit_right,data[int(ind[2]-1):int(ind[4]),1])
#####################################

##### PLOT FIT TO CHECK #####
plt.plot(data[:,1],data[:,2],'kx--')
plt.xlabel('length(Angstrom)')
plt.ylabel('Temperature(K)')
plt.plot(data[int(ind[0]-1):int(ind[2]),1],stretch_left,color='red')
plt.plot(data[int(ind[2]-1):int(ind[4]),1],stretch_right,color='red')
plt.grid(b=None, which='major', axis='both')
#plt.plot(left[:,1],np.polyval(fit_left,left[:,1]),color='red')
#plt.plot(left[:,1],np.polyval(fit_left,left[:,1]),color='red')
plt.show()
##########################

##### WRITE INFO TO A TEXT FILE #####
print('\nFind the temperature fit data in "gradient.txt!"\n')
with open('gradient.txt','w') as fid:
    fid.write('THIS FILE CONTAINS TEMP GRADIENT DATA FOR THE INPUT TEMP PROFILE '
              + filename + '.\n\n')
    fid.write('\tLEFT CHUNK:\n\tslope:\t\t'+str(fit_left[0])+' K/Ang\n\tintercept:\t'
              +str(fit_left[1])+' K\n\n')
    fid.write('\tRIGHT CHUNK:\n\tslope:\t\t'+str(fit_right[0])+' K/Ang\n\tintercept:\t'
              +str(fit_right[1])+' K\n\n')
    fid.write('\tDELTA T:\t'+str(abs(stretch_left[-1]-stretch_right[0]))+' K\n')
        
        
        
        
        
        
        
        
        
        
        

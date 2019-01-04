#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 25 14:58:06 2018

@author: ty
"""
import numpy as np

velsfile = 'example'

steps = 50000 #2000000 #run time
dt = 0.5e-15
dn = 50 #10 #print frequency
prints = steps/dn #times data is printed
split = 1 #times to split data for averaging
tn = prints/split
win = 0.3

time = np.arange(0,tn)*dt*dn*1e6


with open(velsfile, 'r') as fid:
    idsSi = np.zeros((1,1))
    idsGe = np.zeros((1,1))
    for j in range(3):
        fid.readline()
    N = int(fid.readline().strip().split()[0]) #number of atoms
    for j in range(5):
        fid.readline()
    for j in range(N):
        tmp = fid.readline().strip().split()
        if int(tmp[1]) == 1:
            idsSi = np.append(idsSi,j)
        else:
            idsGe = np.append(idsGe,j)
    idsSi = np.delete(idsSi,0)
    idsGe = np.delete(idsGe,0)
    idsSi = idsSi.astype(int)
    idsGe = idsGe.astype(int)
    fid.seek(0)
  
    for i in range(split):
        vels = np.zeros((tn,N*3))
        for j in range(tn):
            data = np.zeros((N,3))
            for k in range(9):
                fid.readline()
            for k in range(N):
                tmp = fid.readline().strip().split()
                data[k,0] = tmp[2]
                data[k,1] = tmp[3]
                data[k,2] = tmp[4]
            vels[j,:] = np.reshape(data,N*3)
            
## GAUSSIAN SMOOTHING ##
#win = round(win*1e12*2*np.pi/0.04) #I'm not certain how this works ...
#gwin = win-1
#n = np.arange(0,gwin+1)-(gwin/2) 
#n = (2.5*n/(gwin/2.0))
#gauss = np.exp(-(1.0/2.0)*np.multiply(n,n)) #gaussian fx to convolve with
#gauss = gauss/np.sum(gauss) #normalized gaussian

#velsG = np.convolve(vels[:,0],gauss,mode='same')

vacf = np.zeros((1000,3))
vavg = np.mean(vels[:,0:3],axis=1)
for i in range(1000):
    x = 1000-i
    
    tmp = np.zeros((x,1))
    for k in range(x):
        tmp[k] = vels[k,0]*vels[k+i,0] 
    vacf[i,0] = tmp.mean()
    
    tmp = np.zeros((x,1))
    for k in range(x):
        tmp[k] = vels[k,0]*vels[k+i,3] 
    vacf[i,1] = tmp.mean()
    
    tmp = np.zeros((x,1))
    for k in range(x):
        tmp[k] = vavg[k]*vavg[k+i] 
    vacf[i,2] = tmp.mean()
    
    
#vacf[:,0]=vacf[:,0]/vacf[:,0].sum(axis=0)
#vacf[:,1]=vacf[:,1]/vacf[:,1].sum(axis=0)
#vacf[:,2]=vacf[:,2]/vacf[:,2].sum(axis=0)
        
        
        
            
        
    
            
            
            
            
            
            
            
            
            
            
            
            
            
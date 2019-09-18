#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 11 15:40:36 2017

@author: michael
"""
import numpy as np

def readlogfile(logfile):
    fin = open(logfile+".log", 'r')
    lines = fin.readlines()
    lines = lines[len(lines)/5:]
    
    store = {}
    datalen = 0
    for line in lines:
        spl = line.split()
        label = spl[0]
        time = float(spl[1])
        cpucount = int(spl[2])
        datalen = int(spl[3])
        
        arr = store.get(label, [])
        arr.append(time)
        
        store[label] = arr
        
    ret = {}
    for key in store:
        ret[key] = (key, np.mean(store[key]), np.std(store[key]), len(store[key]))        
    return ret,datalen

def readcalculations(infile):
    fin = open(infile+".calculations", 'r')
    spl = fin.readlines()[0].split(",")
    return int(spl[2]), int(spl[3])

def readtruncdatalen(infile):
    fin = open(infile+".length", 'r')
    return int(fin.readlines()[0].strip())
        

    
[500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000,6500,7000,7500,8000]
[2.68946,20.4772,67.2049,158.815,314.941,546.643,871.647,1295.37,1865.02,2531.97,3422.53,4617.03,5781.96,6955.75,9030.97,10748.7]
[1.26579,0.594652,1.26787,2.17737,3.34845,4.80933,6.64205,8.55771,10.9824,13.7973,17.2328,20.5683,25.3889,30.3424,35.2276,41.72]
[4.6105,36.6315,124.708,298.266,583.067,994.478,1570.76,2323.91,3314.18,4544.31,6046.77,7842.74,10064.7,12791.8,15235.3,18328.7]
[0.40122,0.92538,2.42256,4.74363,7.91642,11.7285,16.4717,22.0277,28.6255,36.165,44.7627,53.8028,65.2386,77.7947,92.4882,108.137]

    
        
import matplotlib.pyplot as plt

#x = [8,16,32,48,52,56,60,64,68,72,76,80,96,112]
#y2 = [672,1280,2880,4992,5600,6240,6912,7616,8352,9120,9920,10752,14400,18560]
#y1 = [2.21E+06,1944279.893,1857179.299,1807270.051,1803443.1,1796086.342,1794303.73,1775866.719,1776429.371,1778703.051,1782028.475,1790678.554,1836030.025,1847181.208]


plt.rc('text', usetex=True)
plt.rc('font', family='serif')

fig, ax1 = plt.subplots(figsize=(8, 5))
x = [500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000,6500,7000,7500,8000]

inside_cpu_data = [2.68946,20.4772,67.2049,158.815,314.941,546.643,871.647,1295.37,1865.02,2531.97,3422.53,4617.03,5781.96,6955.75,9030.97,10748.7]
inside_cuda_data = [0.301061,0.594652,1.26787,2.17737,3.34845,4.80933,6.64205,8.55771,10.9824,13.7973,17.2328,20.5683,25.3889,30.3424,35.2276,41.72]
outside_cpu_data = [4.6105,36.6315,124.708,298.266,583.067,994.478,1570.76,2323.91,3314.18,4544.31,6046.77,7842.74,10064.7,12791.8,15235.3,18328.7]
outside_cuda_data = [0.40122,0.92538,2.42256,4.74363,7.91642,11.7285,16.4717,22.0277,28.6255,36.165,44.7627,53.8028,65.2386,77.7947,92.4882,108.137]
#inside_cpu_data = outside_cpu_data
#inside_cuda_data = outside_cuda_data

for  (x1, y1,y2) in zip(x,inside_cpu_data,inside_cuda_data):
    print(x1,y1/y2,y2)


ax1.plot(x, inside_cpu_data, 'ro', label="CPU (AMD Ryzen 1700X)")
ax1.plot(x, inside_cuda_data, 'bo', label="GPU (Geforce GTX 1080)")


#ax2.plot([i for i in xrange(1,113)], [12*numHiddenStates + (numHiddenStates-1) + (numHiddenStates*numHiddenStates - numHiddenStates) + numHiddenStates*38 + numHiddenStates*2 for numHiddenStates in xrange(1,113)], 'k')

ymarks = [1.0e0,1.0e1,1.0e2,1.0e3,1.0e4,1.0e5]
ax1.set_xlabel('Number of alignment sites', fontsize=16)
ax1.set_ylabel(r'Mean time in seconds', color='k', fontsize=16)
for  (x1, y1,y2) in zip(x,inside_cpu_data,inside_cuda_data):
    plt.text(x1+60, y2*1.4, r'$%d\times$'% int(y1/y2),fontsize=11,horizontalalignment='center')
#ax2.set_ylabel(r'--- Number of free parameters', color='k', fontsize=16)
plt.title("Inside algorithm timings", fontsize=17)
ax1.tick_params(axis='both', which='major', labelsize=14)
ax1.set_yscale("log")
plt.grid(which='major', axis='both')
ax1.set_yticks(ymarks)
#ax2.tick_params(axis='both', which='major', labelsize=14)
plt.legend()
plt.subplots_adjust(left=0.15, right=0.87, top=0.93, bottom=0.12)
plt.savefig("insidecpugputimeings.svg")
#plt.show()
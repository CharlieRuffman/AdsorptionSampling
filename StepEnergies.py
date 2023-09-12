#!/usr/bin/env python3

import struct
import ase
from ase import io
import numpy as np
import multiprocessing as mp 
import os
import sys

selectEvery=80 #Resolution of sampling in timesteps - must match the value set in plotShadowMD.py

#print what selectEvery is set to
print("Resolution set to take every:", selectEvery, "images (can be changed in the script)")

try:
    filename=sys.argv[1]
except:
    print("No input file name given in first argument")
    exit()


print("Reading input - may take some time if large files")
try:
    a=io.read(filename,":", format="vasp-out") #read in the outcar file
except:
    print("Input file is not a regular vasp OUTCAR, attempting to read as another format. May cause errors with energy values")
    try:
        a=io.read(filename,":")
    except:
        print("Could not read input file - stopping")
        exit()

print("Input file loaded")
a=a[::selectEvery] #select every nth image


io.write("SelectedImg.traj",a)
data=np.zeros((len(a), 2))
data[:,0]=np.arange(0,len(a))
data[:,1]=[b.get_potential_energy() for b in a]


np.savetxt("StepEnergies.txt", data)

print("Done")




        
        
        



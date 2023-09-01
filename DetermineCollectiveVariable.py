#!/usr/bin/env python3
#Python script to take two structures and calculate the RMSD between them using the Kabsch algorithm

import numpy as np
import sys
import ase
from ase import io
from ase import Atoms

#read in the two structures
try:
    struct1=io.read(sys.argv[1])
except:
    print("Could not read first structure")
    exit()

try:
    struct2=io.read(sys.argv[2])
except:
    print("Could not read second structure")
    exit()

print("Files Loaded")
print("Important: ensure atom order is the same in both files")

#check that the structures have the same number of atoms
if len(struct1)!=len(struct2):
    print("Structures have different number of atoms")
    exit()

#calculate the RMSD for each atom and each coordinate separately
xrmsd=np.sqrt(((struct1.get_positions()[:,0]-struct2.get_positions()[:,0])**2))
yrmsd=np.sqrt(((struct1.get_positions()[:,1]-struct2.get_positions()[:,1])**2))
zrmsd=np.sqrt(((struct1.get_positions()[:,2]-struct2.get_positions()[:,2])**2))


#put x y and z RMSD into a single array
rmsd=np.zeros(len(xrmsd)*3)
rmsd[0::3]=xrmsd

#Write the RMSD to a file for each atom including it's element and index
with open("RMSD.txt","w") as f:
    for i in range(len(rmsd)):
        f.write(str(i)+" "+str(struct1[i].symbol)+" "+str(rmsd[i])+"\n")

print("RMSD calculated and written to RMSD.txt")
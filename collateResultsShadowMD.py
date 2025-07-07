#!/usr/bin/env python3


import ase
from ase import io
import numpy as np
import logging
import multiprocessing as mp 
import os
import sys

'''
User set variable for the energy of the adsorbate in a box - has to be calculated separately 
'''
adsorbateEBox=1/2*-6.40408151
#-21.97165359+1/2*-6.40408151 #formate
#-30.205455 - -3.38362605*2 # methanol without H2
#-30.205455 # methanol
#-6.40408151 #H2

#range of image files created by setupShadowMD.py to process
minImage=sys.argv[1] #first image to process
maxImage=sys.argv[2] #last image to process


class ShadowMD():
    def __init__(self,maxImage, minImage=0):
        print("Initialising")
        self.minImage=int(minImage)
        self.maxImage=int(maxImage)


    def run(self):
        #takes a SINGLE structure and performs a single-point energy evaluation
        computedTrajectory=ase.io.Trajectory("OutputTrajectory.traj", mode="a")
        stepEnergies=None
        print("Processing images from "+str(self.minImage)+" to "+str(self.maxImage))
        for i in range(self.minImage,self.maxImage+1):
            try:
                struct=io.read("image"+str(i)+"/OUTCAR", ":")
            except:
                print("structure "+ str(i)+ " not found, continuing")
                continue
            energy=struct[0].get_potential_energy()
            computedTrajectory.write(struct[0])
            if i==self.minImage:
                stepEnergies=np.array([[i, energy, energy+adsorbateEBox]])
            else:
                stepEnergies=np.append(stepEnergies,[[i, energy, energy+adsorbateEBox]], axis=0)
        np.savetxt("StepEnergy.txt",stepEnergies)
     



if __name__ == '__main__':
    shadowMD=ShadowMD(minImage=minImage, maxImage=maxImage)
    shadowMD.run()



        
        
        



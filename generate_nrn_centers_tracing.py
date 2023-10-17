"""
Lauren R. Madden 2020
"""
#########################################################################################################
# Finding compartment centers for COMSOL#################################################################
#########################################################################################################
import os
import numpy as np
from neuron import h
from NEURON.Cell import Ad_Cell

print("Creating model coordinate file for COMSOL...")

# Imported in mm

points_iseg = np.loadtxt('NEURON/XYZcoords/iSeg.txt', delimiter="\t")
points_stem = np.loadtxt('NEURON/XYZcoords/lowerStem.txt', delimiter="\t")
points_per = np.loadtxt('NEURON/XYZcoords/peripheralAxon.txt', delimiter="\t")
points_cen = np.loadtxt('NEURON/XYZcoords/centralAxon.txt', delimiter="\t")

isegDiam = 1.2 # Amir & Devor 2003 Table 1: 5/8 = x/1.9 (MRG fiberD 5.7 nodeD) -> x = 1.2
# Following diams from Lee Chung Chung Coggeshall 1986 J of Comparative Neurology
stemDiam = 5.7
perDiam = 5.7
cenDiam = 3.0

cellAd = Ad_Cell(axon_trajectory_iseg = points_iseg, fiberD_iseg = isegDiam, \
                axon_trajectory_stem = points_stem, fiberD_stem = stemDiam, \
                axon_trajectory_per  = points_per,  fiberD_per  = perDiam, \
                axon_trajectory_cen  = points_cen,  fiberD_cen  = cenDiam)

# Export compartment centers to text file
centers = np.zeros([np.shape(cellAd.get_secs())[0],3])
i = 0
for x in cellAd.get_secs():
   centers[i,:] = np.array(cellAd.center_3Dcoords(x))*1e-3 # convert from um to mm
   i += 1
   
np.savetxt('InterpolationFiles/nrn_centers.txt',centers,delimiter="\t")
"""
Lauren R. Madden 2020
"""
#########################################################################################################
###### Finding compartment centers for COMSOL with varied GIS trajectory ################################
#########################################################################################################

# st = "sp0p1" # spiral and tortuosity
st = "rw0p1" # random walk and tortuosity
# st = "len130" # spiral and length

import os
import sys
import numpy as np
from neuron import h
from NEURON.Cell import Ad_Cell

# Imported in mm
points_stem = np.loadtxt('NEURON/XYZcoords/lowerStem.txt', delimiter="\t")
points_per = np.loadtxt('NEURON/XYZcoords/peripheralAxon.txt', delimiter="\t")
points_cen = np.loadtxt('NEURON/XYZcoords/centralAxon.txt', delimiter="\t")

stem_top_pt = points_stem[0]
points_iseg = np.loadtxt('NEURON/XYZcoords/iSeg_0_'+st+'.txt', delimiter="\t")

# reverse spiral points. original spiral is built bottom to top. iseg needs to be top to bottom.
rev_iseg = np.zeros(np.shape(points_iseg))
for i in range(np.shape(points_iseg)[0]):
    rev_iseg[i] = points_iseg[-(i+1)]
points_iseg = rev_iseg*1e-3 # um to mm

# add points_iseg on top of points_stem
points_iseg[-1,:] = stem_top_pt # make last (bottom) iseg point same as stem top
for i in range(np.shape(points_iseg)[0]-1):
 	points_iseg[i][0] = points_iseg[i][0] + stem_top_pt[0]
 	points_iseg[i][1] = points_iseg[i][1] + stem_top_pt[1]
 	points_iseg[i][2] = points_iseg[i][2] + stem_top_pt[2]
 	
np.savetxt('NEURON/XYZcoords/iSeg_'+st+'.txt',points_iseg,delimiter="\t")

isegDiam = 1.2 # Amir & Devor 2003 Table 1: 5/8 = x/1.9 (MRG fiberD 5.7 nodeD) -> x = 1.2
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
    

np.savetxt('InterpolationFiles/nrn_centers_'+st+'.txt',centers,delimiter="\t")

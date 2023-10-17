"""
Lauren R. Madden, 2023

Glomerulous Trajectory - Spiral with same height

Ends on point directly above starting point at specified distance away
"""
totalHeight = 32

isegDiam = 1.2 # diameter of compartment
n = 100 # number of compartments & length in um
c_len = 1 # compartment length in um
cmpt_rad = isegDiam/2

"""
    Determine the coordinates of the center of glomerulus compartments determined by a spiral algorithm.

    Parameters:
        - totalheight = total global height of the trajectory
        - isegDiam = diameter of the compartment
        - n = number of compartments. each compartment
        - cmpt_len = length of each compartment in um
        - cmpt_rad = radius of each compartment


    Outputs:
        - array of x,y,z coordinates (um) of compartments
    
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

start = [0,0,0]

# build trajectory
pos = np.zeros((n,3)) # position of center of each compartment

r = 11.5 # radius of spiral (11.5 correlates with t = 0.5)
num_ext = 12 # number of compartments extending towards center (12 correlates with t = 0.5)
ext_height = np.sqrt(-(r**2 - num_ext**2)) # height of extension towards center of spiral
h = totalHeight-(2*ext_height)-2 # spiral height

b = np.sqrt((n*c_len)**2 + h**2)/(2*np.pi*r) # number of full rotations (total length of trajectory)
b = b*2*np.pi # number of rotations in radians

# first point is at origin, second point is straight up
pos[1,0] = 0
pos[1,1] = 1
pos[1,2] = 0

theta = np.arcsin(ext_height/num_ext)
for j in range(2,int(num_ext+2)):
    pos[j,0] = pos[j-1,0]+c_len*np.cos(np.pi)*np.sin(theta-(np.pi/2))
    pos[j,1] = pos[j-1,1]-c_len*np.cos(theta+(np.pi/2))
    pos[j,2] = pos[j-1,2]+c_len*np.sin(np.pi)*np.sin(theta-(np.pi/2))

k = j

for i in range(2*int(num_ext+2),2*n-2*int(num_ext+1),2):
    j = int(i/2)
    loc = b*(i-2*int(num_ext+2)+2)/(2*n-2*int(num_ext+1)-2*int(num_ext+2)) # angle in radians of ith compartment
    pos[j,0] = r*np.cos(loc) # x pos
    pos[j,1] = (ext_height+1)+h*(i-2*int(num_ext+2)+2)/(2*n-2*int(num_ext+1)-2*int(num_ext+2)) # y pos
    pos[j,2] = r*np.sin(loc) # z pos

theta = np.arcsin(ext_height/num_ext)
for i in range(2*n-2*int(num_ext+1),2*n-2,2):
    j = int(i/2)
    pos[j,0] = pos[j-1,0]-c_len*np.cos(loc)*np.sin(theta+(np.pi/2))
    pos[j,1] = pos[j-1,1]+c_len*np.cos(theta-(np.pi/2))
    pos[j,2] = pos[j-1,2]-c_len*np.sin(loc)*np.sin(theta+(np.pi/2))


pos[n-1,0] = 0
pos[n-1,1] = 32
pos[n-1,2] = 0

# viewing trajectory
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.plot(pos[:,0],pos[:,2],pos[:,1])
ax.set_xlabel("X pos (um)")
ax.set_ylabel("Z pos (um)")
ax.set_zlabel("Y pos (um)")
ax.set_title("Num cmpts = "+str(n)+", # spirals = "+str(round(b/(2*np.pi),2)))
plt.show()

np.savetxt("NEURON/XYZcoords/iSeg_0_len"+str(n)+".txt",pos,delimiter='\t')
"""
Lauren R. Madden, 2023

Glomerulous Trajectory - Random Walk

Ends on point directly above starting point at specified distance away

"""

t = 0.1 # tortuosity level (between 0 and 1)
totalheight = 40 # total height of the space the path fills

n = 200 # number of compartments & length in um
c_len = 1 # compartment length in um
c_rad = 0.6 # compartment diameter in um

"""
    Determine the coordinates of the center of glomerulus compartments determined by a random walk algorithm.

    Parameters:
        - tort = tortuosity metric, must be between 0 and 1 (least to most tortuous)
        - totalHeight = total global height of the trajectory
        - num = number of compartments. each compartment
        - cmpt_len = length of each compartment in um
        - cmpt_rad = radius of each compartment


    Outputs:
        - array of x,y,z coordinates (um) of compartments
    
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm


if t < 0 or t > 1:
    raise TypeError('Please provide a valid tortuosity value between 0 and 1.')

# measure of tortuosity (0.1 = least tortuous, 0.9 = most tortuous)

pos = np.zeros((n,3)) # position of center of each compartment
taken = np.zeros((1,3)) # list for taken positions

r = -15*t + 21.5 # limit outer radius of the path (max = 20, min = 8)

degrange = (100*t)+50
degmin = 90 - degrange/2
degmax = 90 + degrange/2

if np.mod(10*t,2) != 0:
    num_ext = int((totalheight-2*r)/2) # number of extending compartments per branch
else:
    num_ext = int((totalheight-2*r)/2) + 1


####################
# Helper functions #
####################

# create function for checking if point is within sphere constraint
# x,y,z = coordinates, r = radius of sphere
def iswithinsphere(x,y,z,r):
    sphcen = [0,r+num_ext,0] # point at center of sphere, on y axis above xz origin
    dist = np.sqrt( (sphcen[0]-x)**2 + (sphcen[1]-y)**2 + (sphcen[2]-z)**2) # distance between point and center of sphere
    # return whether point is within radius distance of center
    return dist <= r, dist

# create function for checking if point is outside taken points
# x,y,z = coordinates, taken = taken list, c_rad = compartment diameter
def isoutsidetaken(x,y,z,taken,c_rad):
    if np.shape(taken)[0] > (-3.5*t+3.35)+2:
        for i in range(np.shape(taken)[0]-int(-3.5*t+3.35)-2): # iterate through each point in taken, other than most recent
            dist = np.sqrt( (taken[i,0]-x)**2 + (taken[i,1]-y)**2 + (taken[i,2]-z)**2 ) # get distance between point and taken entry
            if dist <=  (-3.5*t+3.35) + np.sqrt(c_rad**2 + (c_len/2)**2): # if point is too close to a taken entry (less than c_rad plus 1-tort metric), return false
                return False
    return True # if point makes it through taken list, return true

# create a function that returns which point is closest to where the soma should start
# takes in coordinate arrays
def pickclosest(xx,yy,zz):
    somaloc = [0,2*r-c_len+num_ext,0] # point where soma should start
    dist = []
    for j in range(np.shape(xx)[0]):
        dist.append(np.sqrt((xx[j]-somaloc[0])**2 + (yy[j]-somaloc[1])**2 + (zz[j]-somaloc[2])**2))
    return np.where(dist==min(dist))[0][0]

# create a function that returns distance between input point and top of the sphere
# takes in coordinate arrays
def distfromend(x,y,z):
    somaloc = [0,2*r-c_len+num_ext,0] # point where soma should start, top of sphere
    return np.sqrt((x-somaloc[0])**2 + (y-somaloc[1])**2 + (z-somaloc[2])**2)

# create function that converts degrees to radians
# takes in angle in degrees
def rad(deg):
    return deg*np.pi/180


####################
# Build Trajectory #
####################

i = 1

# fill in first several compartments for lower extension
while i <= num_ext:
    pos[i,0] = 0
    pos[i,1] = i
    pos[i,2] = 0
    i+=1

i_hold = i
pos_hold = pos.copy()

ang1_prev = 0
ang2_prev = 0
b = 0

flag_traj = True # flag for trajectories that don't end close enough
maxtries1 = 0
maxtries2 = 0


flag_while = False # flag for stuck while loops

while flag_traj:

    while i < n-num_ext:
        
        if i > num_ext+1:
            nc = i/n # metric for current compartment number along total path
            ed = distfromend(x,y,z)/(2*r) # metric for current distance from top of soma
            b = int(1-(0.2*ed)/(nc-1-0.2*t)) # bias number. number of options to select from to choose point closest to soma start point
            if b < 1: # if bias number is less than one, default 1
                b = 1
            
        if b > 1: # if adding bias
            # initialize empty arrays for selecting from for ith compartment
            xx = np.zeros(b)
            yy = np.zeros(b)
            zz = np.zeros(b)
            for j in range(b): # generate b possible new points
                ang1 = rad(np.random.randint(degmin,degmax)) # generate random angle 1 between 0 and 180
                ang2 = rad(np.random.randint(degmin,degmax)) # generate random angle 2 between 0 and 180
                ang1 = ang1_prev - rad(90) + ang1
                ang2 = ang2_prev - rad(90) + ang2
                # create new point with angles
                xx[j] = pos[i-1,0] + c_len*np.cos(ang1)*np.sin(ang2)
                yy[j] = pos[i-1,1] + c_len*np.sin(ang1)*np.sin(ang2)
                zz[j] = pos[i-1,2] + c_len*np.cos(ang2)
            
                # check if new point is within sphere and outside taken trajectory. if not, generate new point.
                while not(iswithinsphere(xx[j],yy[j],zz[j],r)[0]) or not(isoutsidetaken(xx[j],yy[j],zz[j],taken,c_rad)):
                    degtemp = 0
                    if iswithinsphere(xx[j],yy[j],zz[j],r)[1] > r*(0.1*t+.8): # if new point is appraoching sphere edge
                        degtemp = degrange*0.5 # increase degree range by 50%
                    ang1 = rad(np.random.randint(degmin-degtemp,degmax+degtemp)) # generate random angle
                    ang2 = rad(np.random.randint(degmin-degtemp,degmax+degtemp)) # generate random angle
                    ang1 = ang1_prev - rad(90) + ang1
                    ang2 = ang2_prev - rad(90) + ang2
                    # create new point with angles
                    xx[j] = pos[i-1,0] + c_len*np.cos(ang1)*np.sin(ang2)
                    yy[j] = pos[i-1,1] + c_len*np.sin(ang1)*np.sin(ang2)
                    zz[j] = pos[i-1,2] + c_len*np.cos(ang2)
                    maxtries1+=1
                    if maxtries1 > 1000: # if trajectory gets stuck, break out of loop
                        flag_while = True
                        break
            
            if not(flag_while):
                # select point closest to soma start point
                x = xx[pickclosest(xx,yy,zz)]
                y = yy[pickclosest(xx,yy,zz)]
                z = zz[pickclosest(xx,yy,zz)]
            
        else:
            ang1 = rad(np.random.randint(degmin,degmax)) # generate random angle 1 between 0 and 180
            ang2 = rad(np.random.randint(degmin,degmax)) # generate random angle 2 between 0 and 180
            if i > 1: # if after first compartment, adjust new angle to be relative to previous angle
                ang1 = ang1_prev - rad(90) + ang1
                ang2 = ang2_prev - rad(90) + ang2
            # create new point with angles
            x = pos[i-1,0] + c_len*np.cos(ang1)*np.sin(ang2)
            y = pos[i-1,1] + c_len*np.sin(ang1)*np.sin(ang2)
            z = pos[i-1,2] + c_len*np.cos(ang2)
        
            # check if new point is within sphere and outside taken trajectory. if not, generate new point.
            while not(iswithinsphere(x,y,z,r)[0]) or not(isoutsidetaken(x,y,z,taken,c_rad)):
                degtemp = 0
                if iswithinsphere(x,y,z,r)[1] > r*0.8: # if new point is appraoching sphere edge
                    degtemp = degrange*0.5 # increase degree range
                ang1 = rad(np.random.randint(degmin-degtemp,degmax+degtemp)) # generate random angle
                ang2 = rad(np.random.randint(degmin-degtemp,degmax+degtemp)) # generate random angle
                if i > 1: # if after first compartment, adjust new angle to be relative to previous angle
                    ang1 = ang1_prev - rad(90) + ang1
                    ang2 = ang2_prev - rad(90) + ang2
                # create new point with angles
                x = pos[i-1,0] + c_len*np.cos(ang1)*np.sin(ang2)
                y = pos[i-1,1] + c_len*np.sin(ang1)*np.sin(ang2)
                z = pos[i-1,2] + c_len*np.cos(ang2)
                maxtries2+=1
                if maxtries2 > 1000: # if trajectory gets stuck, break out of loop
                    flag_while = True
                    break
        
        # once sufficient point is found, add to pos array
        # but first, check if trajectory got stuck
        if flag_while:
            pos = pos_hold.copy() # reset pos matrix
            taken = np.zeros((1,3)) # reset list for taken positions
            flag_while = False # reset flag tools and compartment count
            maxtries1 = 0
            maxtries2 = 0
            i = i_hold
            b = 0
        
        else:
            ang1_prev = ang1
            ang2_prev = ang2
            pos[i,0] = x
            pos[i,1] = y
            pos[i,2] = z
            taken = np.append(taken,np.array([[x,y,z]]),axis=0)
            i+=1

    # back out to trajectory while loop    
    if distfromend(pos[i-1,0],pos[i-1,1],pos[i-1,2]) < 8:
        flag_traj = False
    else:
        print("Invalid trajectory... regenerating")
        ang1_prev = 0
        ang2_prev = 0
        b = 0
        pos = pos_hold.copy() # reset pos matrix
        taken = np.zeros((1,3)) # reset list for taken positions
        flag_while = False # reset flag tools and compartment count
        maxtries1 = 0
        maxtries2 = 0
        i = i_hold


posog = pos.copy() # copy original pos

# Meet end to sphere top
flag_dist = True
curpoint = [pos[i-1,0],pos[i-1,1],pos[i-1,2]] # original last point
if np.mod(10*t,2) != 0:
    somaloc = [0,2*r-c_len+num_ext,0] # point before where soma should start
else:
    somaloc = [0,2*r-2*c_len+num_ext,0]
c = 2 # reverse indexer
numchecked = c_len # number of compartments checked, also length in um
while flag_dist:
    dist = np.sqrt( (somaloc[0]-curpoint[0])**2 + (somaloc[1]-curpoint[1])**2 + (somaloc[2]-curpoint[2])**2)
    if dist > numchecked-c_len:
        curpoint = [pos[i-c,0],pos[i-c,1],pos[i-c,2]] # iterate to next point
        c+=1
        numchecked += c_len
    else:
        flag_dist = False
    
anchor = [pos[i-c,0],pos[i-c,1],pos[i-c,2]] # last point of original random walk
vector_mag = np.sqrt( (somaloc[0]-anchor[0])**2 + (somaloc[1]-anchor[1])**2 + (somaloc[2]-anchor[2])**2) # find magnitude of vector between anchor and soma start
vector_dir = [(somaloc[0]-anchor[0])/vector_mag, (somaloc[1]-anchor[1])/vector_mag, (somaloc[2]-anchor[2])/vector_mag] # find unit vector direction

# rebuild last segment of trajectory towards soma start
for ii in range(numchecked-1):
    pos[i-(c-ii-1),0] = pos[i-(c-ii),0] + vector_dir[0]
    pos[i-(c-ii-1),1] = pos[i-(c-ii),1] + vector_dir[1]
    pos[i-(c-ii-1),2] = pos[i-(c-ii),2] + vector_dir[2]
    
# add last compartment to top of sphere

if np.mod(10*t,2) != 0:
    pos[i-1,0] = 0
    pos[i-1,1] = num_ext+2*r
    pos[i-1,2] = 0
    while i < n:
        pos[i,0] = 0
        pos[i,1] = num_ext+2*r+i-(n-num_ext)+1
        pos[i,2] = 0
        i+=1
else:
    pos[i-1,0] = 0
    pos[i-1,1] = num_ext+2*r-c_len
    pos[i-1,2] = 0
    while i < n:
        pos[i,0] = 0
        pos[i,1] = num_ext+2*r+i-(n-num_ext)
        pos[i,2] = 0
        i+=1

posnew = pos.copy()

# viewing trajectory
xmin = min(pos[:,0])
xmax = max(pos[:,0])
ymin = min(pos[:,1])
ymax = max(pos[:,1])
zmin = min(pos[:,2])
zmax = max(pos[:,2])
fig = plt.figure()
cmap = cm.get_cmap('jet',n)
ax = fig.add_subplot(projection='3d')
for j in range(n):
    ax.plot(pos[j:j+2,0],pos[j:j+2,2],pos[j:j+2,1],color=cmap(j))
ax.set_xlabel("X pos (um)")
ax.set_ylabel("Z pos (um)")
ax.set_zlabel("Y pos (um)")
ax.set_title("Tortuosity = "+str(t)+", Radius = "+str(r)+": "+str(int(xmax-xmin))+"x"+str(int(ymax-ymin))+"x"+str(int(zmax-zmin)))
plt.xlim([-r,r])
plt.ylim([-r,r])
plt.show()


np.savetxt("NEURON/XYZcoords/iSeg_0_rw"+str(t).replace('.','p')+".txt",pos,delimiter='\t')
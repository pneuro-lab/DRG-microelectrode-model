def find_devor_compart_coords_high_res(fiberD,axon_trajectory,branch):
    """
        Determine the coordinates of Devor 2003 axon nodes, mysa, flut, and stin for a given high-resolution trajectory.
        Spacing is determined by the fiber diameter, axon branch type, and ratios in Amir & Devor 2003

        Lauren R. Madden 2021

        Inputs:
            - fiberD = um, fiber diameter
            - axon_trajectory = mm, n x 3 array containing the line segments describing the axon trajectory
            - branch = either 'iseg' (initial segment), 'stem_with_iseg', (myelinated stem axon), 'peripheral', 'central'

        Outputs:
            - array of x,y,z coordinates (um) of compartments with the proper spacing along the provided trajectory
            - gives one set of coordinates for each compartment separation

        Assumptions:
            - node spacing is greater than the distance between points in axon_trajectory
            - builds top (near cell body) to bottom for stem axon and medial (T-junction) to lateral (away from
                T-junction) for peripheral/central branches
        
    """

    from numpy import array,sqrt,shape

    axon_trajectory = 1e3*array(axon_trajectory) # convert mm to um

    # specify deltax (i.e. node spacing) based on fiber diameter
    if (fiberD==2.0):
        deltax=117
        flut=10
    if (fiberD==3.0):
        deltax=309
        flut=21
    if (fiberD==5.7):
        deltax=500
        flut=35
    if (fiberD==7.3):
        deltax=750
        flut=38
    if (fiberD==8.7):
        deltax=1000
        flut=40
    if (fiberD==10.0):
        deltax=1150
        flut=46
    if (fiberD==11.5):
        deltax=1250
        flut=50
    if (fiberD==12.8):
        deltax=1350
        flut=54 
    if (fiberD==14.0):
        deltax=1400
        flut=56
    if (fiberD==15.0):
        deltax=1450
        flut=58
    if (fiberD==16.0):
        deltax=1500 
        flut=60
    
    nodelen = 1

    if branch == 'stem_with_iseg':
        
        # determine length of axon
        length = 0
        for p in range(shape(axon_trajectory)[0]-1):
            length = length + sqrt(sum((axon_trajectory[p+1,:]-axon_trajectory[p,:])**2))
        
        # determine internode lengths
        base_inter_lengths = [86,131,169,202] # Devor 2003, added 1,1,1,1 (4 total) to account for node lengths
        total_inter_length = sum(base_inter_lengths)
        percentages = list()
        for n in base_inter_lengths:
            percentages.append(n/total_inter_length) # spacing is proportional to entire stem axon length
        dx_var = [n*(length-2) for n in percentages] # leave space for final node
        
        # determine flut lengths
        flut_var = list()
        for n in dx_var:
            flut_var.append(n*(flut/deltax)) # flut length is proportional to original stem axon internode values
        
        ii,jj = 0,0                         #counters for stepping through points of axon trajectories
        xx,yy,zz = list(),list(),list()     #"interpolated" node locations for NEURON
        P0 = axon_trajectory[ii,:]              #initialize point, P0
        P1 = axon_trajectory[ii+1,:]            #initialize point, P1
        xx.append(P0[0])
        yy.append(P0[1])
        zz.append(P0[2]) #define initial point of axon, i.e. first node
        temp = 0
        compt_order = ['node', 'mysa1', 'flut1', 'stin1', 'stin2', 'stin3', 'stin4', 'stin5', 'stin6', 'flut2', 'mysa2']
        flag = 'node' # flag for compartment type
        nodecounter = 0 # counter number of Devor internodes passed
        while True:
            if nodecounter < 5:
                if nodecounter < 4:
                    dx = dx_var[nodecounter] # set total length between nodes
                    flutlen = flut_var[nodecounter] # set flut lengths for stem
                    mysalen = 1 # set mysa length for stem
                    stinlen = (dx - nodelen - (2*mysalen) - (2*flutlen))/6 # set stin lengths for stem
                    
                    # create dictionary with compartment flags and lengths
                    compt = {'node':nodelen,\
                             'mysa1':mysalen,\
                             'flut1':flutlen,\
                             'stin1':stinlen,\
                             'stin2':stinlen,\
                             'stin3':stinlen,\
                             'stin4':stinlen,\
                             'stin5':stinlen,\
                             'stin6':stinlen,\
                             'flut2':flutlen,\
                             'mysa2':mysalen }
                
            else:
                raise TypeError('Check stem length calculations')
            temp = temp + sqrt(sum((P1-P0)**2))
            if (temp < compt[flag]): # if temp still shorter than compt length, check next point
                ii+=1
                if (ii == (shape(axon_trajectory)[0]-1)): # UNLESS: if next pointer is end of axon trajectory
                    if abs((compt[flag] - temp)/compt[flag]) < 0.001: # if end of axon trajectory is close enough to next compartment point
                        xx.append(axon_trajectory[-1,0]) # add last point to output
                        yy.append(axon_trajectory[-1,1])
                        zz.append(axon_trajectory[-1,2])
                    break
                P0 = P1
                P1 = axon_trajectory[ii+1,:]
            elif (temp == compt[flag]) : # if temp equals compartment length, add point to output
                xx.append(P1[0])
                yy.append(P1[1])
                zz.append(P1[2])
                ii+=1
                jj+=1
                if (ii == (shape(axon_trajectory)[0]-1)):
                    break
                P0 = [xx[jj],yy[jj],zz[jj]]
                P1 = axon_trajectory[ii+1,:]
                temp = 0
                if compt_order.index(flag) < 10: # change compartment
                    flag = compt_order[compt_order.index(flag)+1]
                else:
                    flag = 'node'
                    nodecounter = nodecounter+1
            else: #if next point exceeds compartment length, add interpolated point to output
                P0P1 = sqrt(sum((P1-P0)**2))
                tt = (temp-compt[flag])/P0P1
                xx.append((tt)*P0[0]+(1-tt)*P1[0])
                yy.append((tt)*P0[1]+(1-tt)*P1[1])
                zz.append((tt)*P0[2]+(1-tt)*P1[2])
                jj+=1
                ii+=1
                if (ii == (shape(axon_trajectory)[0]-1)):
                    break
                P0 = [xx[jj],yy[jj],zz[jj]]
                P1 = axon_trajectory[ii+1,:]
                temp = 0
                if compt_order.index(flag) < 10: # change compartment
                    flag = compt_order[compt_order.index(flag)+1]
                else:
                    flag = 'node'
                    nodecounter = nodecounter+1
        return array([xx,yy,zz]).T
    
    elif branch == 'peripheral':
        
        dx_var = [x*deltax for x in [461/1567, 670/1567, 1119/1567]] # Amir & Devor 2003
        
        # determine flut lengths
        flut_var = list()
        for n in dx_var:
            flut_var.append(n*(flut/deltax)) # flut length is proportional to original stem axon internode values
        
        ii,jj = 0,0                         #counters for stepping through points of axon trajectories
        xx,yy,zz = list(),list(),list()     #"interpolated" node locations for NEURON
        P0 = axon_trajectory[ii,:]              #initialize point, P0
        P1 = axon_trajectory[ii+1,:]            #initialize point, P1
        xx.append(P0[0])
        yy.append(P0[1])
        zz.append(P0[2]) #define initial point of axon, i.e. first node
        temp = 0
        compt_order = ['node', 'mysa1', 'flut1', 'stin1', 'stin2', 'stin3', 'stin4', 'stin5', 'stin6', 'flut2', 'mysa2']
        flag = 'node' # flag for compartment type
        nodecounter = 0
        while True:
            if nodecounter < 3: # first three internodes are Amir&Devor lengths
                dx = dx_var[nodecounter] # set total length between nodes
                flutlen = flut_var[nodecounter] # set flut lengths for stem
                mysalen = 3 # set mysa length for stem
                stinlen = (dx - nodelen - (2*mysalen) - (2*flutlen))/6 # set stin lengths for stem
                
                # create dictionary with compartment flags and lengths
                # create dictionary with compartment flags and lengths
                compt = {'node':nodelen,\
                         'mysa1':mysalen,\
                         'flut1':flutlen,\
                         'stin1':stinlen,\
                         'stin2':stinlen,\
                         'stin3':stinlen,\
                         'stin4':stinlen,\
                         'stin5':stinlen,\
                         'stin6':stinlen,\
                         'flut2':flutlen,\
                         'mysa2':mysalen }
                
            else: # all other internodes
                dx = deltax # set total length between nodes
                flutlen = flut # set flut lengths for stem
                mysalen = 3 # set mysa length for stem
                stinlen = (dx - nodelen - (2*mysalen) - (2*flutlen))/6 # set stin lengths for stem
                
                # create dictionary with compartment flags and lengths
                # create dictionary with compartment flags and lengths
                compt = {'node':nodelen,\
                         'mysa1':mysalen,\
                         'flut1':flutlen,\
                         'stin1':stinlen,\
                         'stin2':stinlen,\
                         'stin3':stinlen,\
                         'stin4':stinlen,\
                         'stin5':stinlen,\
                         'stin6':stinlen,\
                         'flut2':flutlen,\
                         'mysa2':mysalen }
    
            temp = temp + sqrt(sum((P1-P0)**2))
            if (temp < compt[flag]): # if temp still shorter than compt length, go to next point
                ii+=1
                if (ii == (shape(axon_trajectory)[0]-1)): # UNLESS: if next pointer is end of axon trajectory
                    if abs((compt[flag] - temp)/compt[flag]) < 0.001: # if end of axon trajectory is close enough to next compartment point
                        xx.append(axon_trajectory[-1,0]) # add last point to output
                        yy.append(axon_trajectory[-1,1])
                        zz.append(axon_trajectory[-1,2])
                    break
                P0 = P1
                P1 = axon_trajectory[ii+1,:]
            elif (temp == compt[flag]) : # if temp equals compartment length, add point to output
                xx.append(P1[0])
                yy.append(P1[1])
                zz.append(P1[2])
                ii+=1
                jj+=1
                if (ii == (shape(axon_trajectory)[0]-1)):
                    break
                P0 = [xx[jj],yy[jj],zz[jj]]
                P1 = axon_trajectory[ii+1,:]
                temp = 0
                if compt_order.index(flag) < 10: # change compartment
                    flag = compt_order[compt_order.index(flag)+1]
                else:
                    flag = 'node'
                    nodecounter = nodecounter+1
            else: #if next point exceeds compartment length, add interpolated point to output
                P0P1 = sqrt(sum((P1-P0)**2))
                tt = (temp-compt[flag])/P0P1
                xx.append((tt)*P0[0]+(1-tt)*P1[0])
                yy.append((tt)*P0[1]+(1-tt)*P1[1])
                zz.append((tt)*P0[2]+(1-tt)*P1[2])
                jj+=1
                ii+=1
                if (ii == (shape(axon_trajectory)[0]-1)):
                    break
                P0 = [xx[jj],yy[jj],zz[jj]]
                P1 = axon_trajectory[ii+1,:]
                temp = 0
                if compt_order.index(flag) < 10: # change compartment
                    flag = compt_order[compt_order.index(flag)+1]
                else:
                    flag = 'node'
                    nodecounter = nodecounter+1
        return array([xx,yy,zz]).T
    
    elif branch == 'central':
        
        dx_var = [x*deltax for x in [358/1450, 780/1450, 1170/1450]] # Devor 2003
        
        # determine flut lengths
        flut_var = list()
        for n in dx_var:
            flut_var.append(n*(flut/deltax)) # flut length is proportional to original stem axon internode values
        
        ii,jj = 0,0                         #counters for stepping through points of axon trajectories
        xx,yy,zz = list(),list(),list()     #"interpolated" node locations for NEURON
        P0 = axon_trajectory[ii,:]              #initialize point, P0
        P1 = axon_trajectory[ii+1,:]            #initialize point, P1
        xx.append(P0[0])
        yy.append(P0[1])
        zz.append(P0[2]) #define initial point of axon, i.e. first node
        temp = 0
        compt_order = ['node', 'mysa1', 'flut1', 'stin1', 'stin2', 'stin3', 'stin4', 'stin5', 'stin6', 'flut2', 'mysa2']
        flag = 'node' # flag for compartment type
        nodecounter = 0
        while True:
            if nodecounter < 3: # first three internodes are Amir&Devor lengths
                dx = dx_var[nodecounter] # set total length between nodes
                flutlen = flut_var[nodecounter] # set flut lengths for stem
                mysalen = 3 # set mysa length for stem
                stinlen = (dx - nodelen - (2*mysalen) - (2*flutlen))/6 # set stin lengths for stem
                
                # create dictionary with compartment flags and lengths
                compt = {'node':nodelen,\
                         'mysa1':mysalen,\
                         'flut1':flutlen,\
                         'stin1':stinlen,\
                         'stin2':stinlen,\
                         'stin3':stinlen,\
                         'stin4':stinlen,\
                         'stin5':stinlen,\
                         'stin6':stinlen,\
                         'flut2':flutlen,\
                         'mysa2':mysalen }
                
            else: # all other internodes
                dx = deltax # set total length between nodes
                flutlen = flut # set flut lengths for stem
                mysalen = 3 # set mysa length for stem
                stinlen = (dx - nodelen - (2*mysalen) - (2*flutlen))/6 # set stin lengths for stem
                
                # create dictionary with compartment flags and lengths
                compt = {'node':nodelen,\
                         'mysa1':mysalen,\
                         'flut1':flutlen,\
                         'stin1':stinlen,\
                         'stin2':stinlen,\
                         'stin3':stinlen,\
                         'stin4':stinlen,\
                         'stin5':stinlen,\
                         'stin6':stinlen,\
                         'flut2':flutlen,\
                         'mysa2':mysalen }
    
            temp = temp + sqrt(sum((P1-P0)**2))
            if (temp < compt[flag]): # if temp still shorter than compt length, go to next point
                ii+=1
                if (ii == (shape(axon_trajectory)[0]-1)): # UNLESS: if next pointer is end of axon trajectory
                    if abs((compt[flag] - temp)/compt[flag]) < 0.001: # if end of axon trajectory is close enough to next compartment point
                        xx.append(axon_trajectory[-1,0]) # add last point to output
                        yy.append(axon_trajectory[-1,1])
                        zz.append(axon_trajectory[-1,2])
                    break
                P0 = P1
                P1 = axon_trajectory[ii+1,:]
            elif (temp == compt[flag]) : # if temp equals compartment length, add point to output
                xx.append(P1[0])
                yy.append(P1[1])
                zz.append(P1[2])
                ii+=1
                jj+=1
                if (ii == (shape(axon_trajectory)[0]-1)):
                    break
                P0 = [xx[jj],yy[jj],zz[jj]]
                P1 = axon_trajectory[ii+1,:]
                temp = 0
                if compt_order.index(flag) < 10: # change compartment
                    flag = compt_order[compt_order.index(flag)+1]
                else:
                    flag = 'node'
                    nodecounter = nodecounter+1
            else: #if next point exceeds compartment length, add interpolated point to output
                P0P1 = sqrt(sum((P1-P0)**2))
                tt = (temp-compt[flag])/P0P1
                xx.append((tt)*P0[0]+(1-tt)*P1[0])
                yy.append((tt)*P0[1]+(1-tt)*P1[1])
                zz.append((tt)*P0[2]+(1-tt)*P1[2])
                jj+=1
                ii+=1
                if (ii == (shape(axon_trajectory)[0]-1)):
                    break
                P0 = [xx[jj],yy[jj],zz[jj]]
                P1 = axon_trajectory[ii+1,:]
                temp = 0
                if compt_order.index(flag) < 10: # change compartment
                    flag = compt_order[compt_order.index(flag)+1]
                else:
                    flag = 'node'
                    nodecounter + nodecounter+1
        return array([xx,yy,zz]).T
    
    elif branch == 'iseg':
        
        dx = 1 # node length
        
        ii,jj = 0,0                         #counters for stepping through points of axon trajectories
        xx,yy,zz = list(),list(),list()     #"interpolated" node locations for NEURON
        P0 = axon_trajectory[ii,:]              #initialize point, P0
        P1 = axon_trajectory[ii+1,:]            #initialize point, P1
        xx.append(P0[0])
        yy.append(P0[1])
        zz.append(P0[2]) #define initial point of axon, i.e. first node
        temp = 0
        while True:
            temp = temp + sqrt(sum((P1-P0)**2))
            if (temp < dx):
                ii+=1
                if (ii == (shape(axon_trajectory)[0]-1)):
                    break
                P0 = P1
                P1 = axon_trajectory[ii+1,:]
            elif (temp == dx):
                xx.append(P1[0])
                yy.append(P1[1])
                zz.append(P1[2])
                ii+=1
                jj+=1
                if (ii == (shape(axon_trajectory)[0]-1)):
                    break
                P0 = [xx[jj],yy[jj],zz[jj]]
                P1 = axon_trajectory[ii+1,:]
                temp = 0
            else:
                P0P1 = sqrt(sum((P1-P0)**2))
                tt = (temp-dx)/P0P1
                xx.append((tt)*P0[0]+(1-tt)*P1[0])
                yy.append((tt)*P0[1]+(1-tt)*P1[1])
                zz.append((tt)*P0[2]+(1-tt)*P1[2])
                jj+=1
                ii+=1
                if (ii == (shape(axon_trajectory)[0]-1)):
                    break
                P0 = [xx[jj],yy[jj],zz[jj]]
                P1 = axon_trajectory[ii+1,:]
                temp = 0
        return array([xx,yy,zz]).T
    else:
        raise TypeError('Need to enter correct branch type')
from __future__ import division

if __name__=="__main__":
        
    import sys
    import pickle
    import numpy as np
    import os
    import scipy.io as scio

    PyRoot = os.getcwd().replace('\\','/')

    from NEURON.Cell import Ad_Cell
    from NEURON.rec_dict import rec_dict
    from neuron import h
    import neuron as nrn
    h.load_file("stdrun.hoc")

    import pickle
    import time
    import matplotlib.pyplot as plt

    if str(sys.argv[1]) == 'extra':

        #############################################################################
        ## Define the necessary inputs
        #############################################################################

        interpfile = str(sys.argv[3])
        recDetails = interpfile.split("interpolation")[1]
        recDetails = recDetails.split('.')[0]
        sigmaval = recDetails.split('sigma')[1].split('_')[0].replace('p','.')
        sigmaval = float(sigmaval)

        f_Kmat = "InterpolationFiles/"+str(sys.argv[3]) # interpolation file from COMSOL
        f_save = "resultfile"+recDetails # simulation results pickle file

        #############################################################################
        ## Read in pre-determined K matrix
        #############################################################################
        #read in predetermined K matrix
        #   K matrix has the form [nx ny nz K]
        #   nx,ny,nz = x,y,z coordinates of the center of each compartment
        #   K = the K value for the corresponding compartment
        
        textfile = open(f_Kmat,'r+') 
        output = textfile.readlines()
        output_array = np.empty((np.shape(output)[0],4))
        for i in range(np.shape(output)[0]):
            templist = output[i].split(' ')
            templist = list(filter(lambda a: a != templist[1], templist))
            templist[3] = templist[3][0:-1]
            output_array[i][0] = float(templist[0])*1e3 # convert mm to um
            output_array[i][1] = float(templist[1])*1e3 # convert mm to um
            output_array[i][2] = float(templist[2])*1e3 # convert mm to um
            output_array[i][3] = float(templist[3])
        textfile.close()
        K_mat = output_array

        K = K_mat[:,3]

    elif str(sys.argv[1]) == 'intra':
        GIS = str(sys.argv[2])
        GIS = GIS.split('.')[0]
        f_save = "resultfile_"+GIS # simulation results pickle file
        sigmaval = float(sys.argv[3])
    else:
        raise ValueError("Please specify 'extra' for extracellular voltage waveform or 'intra' for intracellular voltage waveform")


    #############################################################################
    # Define cell model
    #############################################################################
    
    # Imported in mm
    points_iseg = np.loadtxt(PyRoot+'/NEURON/XYZcoords/'+str(sys.argv[2]), delimiter="\t")
    points_stem = np.loadtxt(PyRoot+'/NEURON/XYZcoords/lowerStem.txt', delimiter="\t")
    points_per = np.loadtxt(PyRoot+'/NEURON/XYZcoords/peripheralAxon.txt', delimiter="\t")
    points_cen = np.loadtxt(PyRoot+'/NEURON/XYZcoords/centralAxon.txt', delimiter="\t")
    
    isegDiam = 1.2
    stemDiam = 5.7
    perDiam = 5.7
    cenDiam = 3.0
    
    # load individual axon trajectory
    cell = Ad_Cell(axon_trajectory_iseg = points_iseg, fiberD_iseg = isegDiam, \
                 axon_trajectory_stem = points_stem, fiberD_stem  = stemDiam, \
                 axon_trajectory_per  = points_per,  fiberD_per  = perDiam, \
                 axon_trajectory_cen  = points_cen,  fiberD_cen  = cenDiam)

    #############################################################################
    # define the intracellular stimulus
    #############################################################################

    def test_pulse(amplitude,delay,duration):
        test_pulse = h.IClamp(0.5, sec = cell.get_nodesP()[-2]) # stimulate second to last node of peripheral axon
        test_pulse.amp = amplitude
        test_pulse.delay = delay
        test_pulse.dur = duration
        return test_pulse
                            
    #############################################################################
    # functions to initiliaze and run a NEURON simulation
    #############################################################################

    fs = 500e3
    dt = 1/fs                           # seconds
    tstop = 30                          # ms
    delay = 10                          # ms
    amplitude = 1.0                     # nA
    stim_duration = 0.2                 # ms
    time_duration = tstop               # ms

    h.dt = dt*1e3
    h.celsius = 37.0
    h.tstop = tstop
    h.v_init = -52.0        #mV, initial voltage to begin simulations
    h.steps_per_ms = fs*1e-3


    # --------------------------------------------------
    # Build the Stimulus
    # --------------------------------------------------
    test_pulse = test_pulse(amplitude,delay,stim_duration)
    

    #############################################################################
    ## Create the recording dictionary to keep track of results
    #############################################################################


    recording_dict = dict() # the dictionary that'll be filled with 'keyName':recording pairs
    location = 0.5          # record from the middle of the compartment
    if str(sys.argv[1]) == 'extra':
        variable = 'i_membrane' # the variable to record
    else:
        variable = 'v'

    if variable == 'i_membrane' or variable == 'v':
        for sec in cell.get_secs():
            rec_dict(recording_dict,sec,location,variable)
    else:
        nodesecs = [sec for sec in cell.get_nodesP()]
        nodesecs.extend([sec for sec in cell.get_nodesC()])
        nodesecs.extend([sec for sec in cell.get_nodesT()])
        nodesecs.extend([sec for sec in cell.get_iseg()])
        nodesecs.append([sec for sec in cell.get_soma()])

        for sec in nodesecs:
            rec_dict(recording_dict,sec,location,variable)

    cell.record(recording_dict) # calls the record function from the cell class to actually record the currents


    #############################################################################
    ## Set up Ephaptic Coupling
    #############################################################################

    iseg_sec_to_idx = dict() # dictionary for iseg sections where key:value = section:idx
    iseg_idx = 0
    for sec in cell.get_iseg():
        iseg_sec_to_idx[sec] = iseg_idx
        iseg_idx += 1

    
    # Make N x 3 x 3 Matrix of SectionsIndex x h/center/l x 3D coordinates

    num_iseg = np.shape(cell.get_iseg())[0] # number of iseg sections
    iseg_coords = np.zeros((num_iseg,3,3)) # array for SectionIndex x 3 (h,c,l) x 3 (XYZ)
    for sec in cell.get_iseg():
        iseg_coords[iseg_sec_to_idx[sec]][0][0:3] = cell.retrieve_coordinates(sec)[0][0:3]  # h coords (um)
        iseg_coords[iseg_sec_to_idx[sec]][1][0:3] = list(cell.center_3Dcoords(sec))         # c coords (um)
        iseg_coords[iseg_sec_to_idx[sec]][2][0:3] = cell.retrieve_coordinates(sec)[1][0:3]  # l coords (um)
    

    # Make N x N x 2 Matrix for distances between Generating Section x Receiving Section x h/l
    
    iseg_dists = np.zeros((num_iseg,num_iseg,2)) # array for generating x recieiving x 2 (h or l)
    g_idx = 0 # generating index
    for g_idx in range(num_iseg):
        r_idx = 0 # receiving index
        for r_idx in range(num_iseg):
            iseg_dists[g_idx][r_idx][0] = np.sqrt( (iseg_coords[g_idx][0][0] - iseg_coords[r_idx][1][0])**2 + \
                (iseg_coords[g_idx][0][1] - iseg_coords[r_idx][1][1])**2 + \
                (iseg_coords[g_idx][0][2] - iseg_coords[r_idx][1][2])**2) # calculate gen_rec_h (um)
            iseg_dists[g_idx][r_idx][1] = np.sqrt( (iseg_coords[g_idx][2][0] - iseg_coords[r_idx][1][0])**2 + \
                (iseg_coords[g_idx][2][1] - iseg_coords[r_idx][1][1])**2 + \
                (iseg_coords[g_idx][2][2] - iseg_coords[r_idx][1][2])**2) # calculate gen_rec_l (um)
            r_idx += 1
        g_idx += 1

    sigma_ex = sigmaval # S/m , extracellular conductivity
    rad_iseg = isegDiam/2 # um , iseg radius

    ec_voltages = np.zeros((num_iseg,num_iseg)) # temp matrix for generating x receiving voltages


    ##############################################################################
    ### Run with Ephaptic Coupling
    ##############################################################################

    h.secondorder = 1
    h.load_file("stdrun.hoc")

    print("Simulation Starting...")
    timer = range(0,55,5)
    tt = 0

    while h.t < h.tstop:

        # Calculate ec voltage for each section
        g_idx = 0
        for g_sec in cell.get_iseg():
            r_idx = 0
            for r_sec in cell.get_iseg():
                # get extracellular voltages from generating sections currents and put in temp matrix
                if g_idx != r_idx:
                    ec_voltages[g_idx][r_idx] = ((g_sec.i_membrane * rad_iseg)/(2*sigma_ex)) * np.log(abs( \
                        (np.sqrt(iseg_dists[g_idx][r_idx][0]**2 + rad_iseg**2) - iseg_dists[g_idx][r_idx][0]) / \
                        (np.sqrt(iseg_dists[g_idx][r_idx][1]**2 + rad_iseg**2) - iseg_dists[g_idx][r_idx][1]) )) *\
                        1e-2 # converts ((mA/cm2)*um)/(S/m) = mV 1e-2 -> mV
                else:
                    ec_voltages[g_idx][r_idx] = 0
                r_idx += 1
            g_idx += 1
                
        # Inject ec of all other sections into each section via extracellular mechanism
        for r_sec in cell.get_iseg():
           r_sec.e_extracellular = sum(ec_voltages[:,iseg_sec_to_idx[r_sec]]) # mV

        if (h.t > timer[tt]):
            print("Simulated Time: "+str(timer[tt])+" ms")
            tt+=1

        h.fadvance()

    h.run()

    ###########################################################
    ## Recording Functions
    ###########################################################

    # Downsampling
    def downsample(sig, og_fs, samp_fs):
        samp_bin = range(0,len(sig),int(len(sig)/np.ceil(len(sig)*(samp_fs/og_fs))))
        new_sig = sig[samp_bin]
        return new_sig

    # Theorem of Reciprocity
    def record_recipro():
        # Uses the Theory of Reciprocity to calculate the sfap at each recording electrode specified    
        count=0
        
        summed_compart = np.zeros(len(np.asarray(cell._recordings['t']))) # all compartments summed onto electrode for each timestep

        i_mems = np.zeros((len(cell.get_secs()),len(np.asarray(cell._recordings['t'])))) # individual compartment voltages for each timestep, adjusted by K matrix
        i_mems_cc = np.zeros((len(cell.get_secs()),len(np.asarray(cell._recordings['t'])))) # individual compartment currents for each timestep
        section_rows = dict() # section name dictionary for indexing

        for section in cell.get_secs(): # iterate through each section's recordings
            diam=section.diam # [um]
            dist=section.L # [um]
            circ = np.pi*diam # circumference # [um]
            surfarea = dist*circ # [um^2]
            compartment_current = np.asarray(cell._recordings[section.name()+'(0.5).i_membrane'])*surfarea*1e-8 # this will show the area created to get a measurement of the actual current 
            # currents are saved [mA/cm^2*um^2*10^-8]=[mA] in mA units
            summed_compart = summed_compart + K[count]*compartment_current # add to total waveform
            
            i_mems[count] = K[count]*np.transpose(np.asarray(cell._recordings[section.name()+'(0.5).i_membrane']))*surfarea*1e-8 # individual compartment voltages for each timestep, adjusted by K matrix
            i_mems_cc[count] = np.transpose(np.asarray(cell._recordings[section.name()+'(0.5).i_membrane']))*surfarea*1e-8 # individual compartment currents for each timestep
            
            section_rows[section.name()] = count # record section name index
            count+=1

        return np.asarray(i_mems),section_rows,np.asarray(summed_compart),np.asarray(i_mems_cc)


    def record_voltages():
        # Records compartmental voltages    
        count=0

        v_array = np.zeros((len(cell.get_secs()),len(np.asarray(cell._recordings['t']))))
        section_rows = dict()

        for section in cell.get_secs(): # iterate through each section's recordings
            v_array[count] = np.transpose(np.asarray(cell._recordings[section.name()+'(0.5).v']))
            section_rows[section.name()] = count

            count+=1
        return v_array,section_rows


    def record_ion():
        # Records compartmental ion currents
        count=0

        ion_array = np.zeros((len(nodesecs),len(np.asarray(cell._recordings['t']))))
        section_rows = dict()

        for section in nodesecs: # iterate through each section's recordings
            ion_array[count] = np.transpose(np.asarray(cell._recordings[section.name()+'(0.5).'+variable]))
            section_rows[section.name()] = count

            count+=1
        return ion_array,section_rows   


    # runs appropriate function for recording variable
    if variable == 'i_membrane':
        i_mems,section_rows,waveform,i_mems_cc = record_recipro()
        waveform_ds = downsample(waveform, fs, 30e3)
    elif variable == 'v':
        v_array,section_rows = record_voltages()
    else:
        ion_array,section_rows = record_ion()


    #############################################################################
    # save the desired files
    #############################################################################

    # output the membrane potentials/i_membrane

    # output other values of interest, threshold, AP times and numbers, fiber diameter, stimulus parameters
    if variable == 'i_membrane':
        Pulse_Dict = {'t' : cell._recordings['t'], \
            'waveform' : waveform, \
            'waveform_ds' : waveform_ds, \
            'i_mems' : i_mems, \
            'section_rows' : section_rows}
            # t = ms, waveform = uV, waveform_ds = uV, i_mems = mA/cm^2
    elif variable == 'v':
        Pulse_Dict = {'t' : cell._recordings['t'], \
            'v_array' : v_array, \
            'section_rows' : section_rows}
            # t = ms, v_array = mV
    else:
        Pulse_Dict = {'t' : cell._recordings['t'], \
            'ion_array' : ion_array, \
            'section_rows' : section_rows}
            # t = ms, ion_array = mA/cm^2

    # save results
    if os.path.isdir(PyRoot+'/Results'):
        output = open("Results/"+f_save+"_"+variable,'wb')
    else:
        os.makedirs(PyRoot+'/Results')
        output = open("Results/"+f_save+"_"+variable,'wb')
    pickle.dump(Pulse_Dict,output,protocol=4)
    output.close()


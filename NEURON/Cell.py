from __future__ import division
import os
from neuron import h
import neuron as nrn
h.load_file("stdrun.hoc")
from numpy import pi,shape,array,ones,zeros,sqrt
from NEURON.find_devor_compart_coords_high_res import find_devor_compart_coords_high_res

class Cell(object):
    '''
    A NEURON multi-compartment model
    '''

    def __init__(self,**kwargs):
        """Initialize cell by loading HOC file"""

        self.variables = kwargs #possible variables to define

        #-----------------------------------------------------------------------
        # Cell Parameters
        #-----------------------------------------------------------------------
 
        #Check for directory containing cell model
        if not self.CELL_DIR:
            raise TypeError("No directory specified for loading model axon...")

        #Check for HOC filename of cell model
        if not self.CELL_NAME:
            raise TypeError("No HOC file specified...")

        #-----------------------------------------------------------------------
        # Import Mechanisms
        #-----------------------------------------------------------------------
        nrn.NRN_NMODL_PATH=self.CELL_DIR
        nrn.load_mechanisms(nrn.NRN_NMODL_PATH)
        #----------------------------------------------------------------------
        # Construct cell
        #-----------------------------------------------------------------------
        self._construct_cell()
        #-----------------------------------------------------------------------
        # Load extracellular mechanisms
        #-----------------------------------------------------------------------
        self._xtra()

    def __str__(self):
        return "Generic Axon Type"

    def __del__(self):
        pass

    def _construct_cell(self):
        """ Choose axon type """
        raise "Dummy Axon"


    def center_3Dcoords(self,sec):
        """
        Returns the 3D coordinates of the center of the currently accessed section

        Requires the 'xtra' mechanism
        """
        x = sec(0.5).x_xtra
        y = sec(0.5).y_xtra
        z = sec(0.5).z_xtra

        return x,y,z

    def get_im(self,sec):
        ''' Dummy function for defining im of a section, to be used in 'build_tree' as 'func'
        '''
        im = sec(0.5).i_membrane
        
        return im

    def retrieve_coordinates(self,sec):
        """
        Returns the 3D coordinates of all points in the currently accessed section.

        Does not calculate the center of each section (see 'center_3Dcoords' above)
        """
        xyzds = []
        for ii in range(int(h.n3d(sec=sec))):
            xyzds.append([h.x3d(ii,sec=sec),
                          h.y3d(ii,sec=sec),
                          h.z3d(ii,sec=sec),
                          h.diam3d(ii,sec=sec)])
        return xyzds

    def record(self,recording_dict):
        """ add a new recording for section
            recording_dict should be a dictionary
            e.g.
                recording_dict={'node[0](0.5).v':h.node[0](0.5)._ref_v,
                                'node[20](0.5).v':h.node[20](0.5)._ref_v,
                                'node[40](0.5).v':h.node[40](0.5)._ref_v,
                                'stim.intensity':h.stim._ref_intensity}
        """
        self._recordings = {'t':h.Vector()}
        self._recordings['t'].record(h._ref_t)
        for k,v in recording_dict.items():
            self._recordings[k] = h.Vector()
            self._recordings[k].record(v)

    
    def _xtra(self):
        """ Insert xtra mechanism to use for finding the center of each compartment (see interpxyz.hoc),
            and necessary extracellular mechanisms
        """
        
        for sec in h.allsec():
            sec.insert("xtra")

        #load required xtra files
        h.load_file(self.CELL_DIR+"interpxyz.hoc")    #only interpolates sections that have extracellular
        h.load_file(self.CELL_DIR+"setpointers.hoc")  #automatically calls grindaway() in interpxyz.hoc


    #make setting and getting variables/attributes scalable
    def set_variable(self,k,v):
        self.variables[k] = v
        #k = key, v = value
        
    def get_variable(self,k):
        return self.variables.get(k,None)


        
class Ad_Cell(Cell):
    '''
        A-Delta Fiber

        Lauren R. Madden 2023
    '''
    def __init__(self,**kwargs):

        self.variables = kwargs

        if 'axon_trajectory_iseg' not in self.variables: raise TypeError('Need to specify initial segment trajectory!!!')
        if 'axon_trajectory_stem' not in self.variables: raise TypeError('Need to specify stem axon trajectory!!!')
        if 'axon_trajectory_per' not in self.variables: raise TypeError('Need to specify peripheral axon trajectory!!!')
        if 'axon_trajectory_cen' not in self.variables: raise TypeError('Need to specify central axon trajectory!!!')
        
        if 'fiberD_iseg' not in self.variables:          raise TypeError('Need to specify initial segment diameter!!!')
        if 'fiberD_stem' not in self.variables:          raise TypeError('Need to specify stem axon diameter!!!')
        if 'fiberD_per' not in self.variables:          raise TypeError('Need to specify peripheral axon diameter!!!')
        if 'fiberD_cen' not in self.variables:          raise TypeError('Need to specify central axon diameter!!!')
        
        self.CELL_DIR = os.getcwd().replace('\\','/') + '/NEURON/ADelta_files/'
        self.CELL_NAME = "ADeltaFiber.hoc"
        
        super(Ad_Cell,self).__init__(**kwargs)
       
    def __str__(self):
        return "A-Delta Fiber"


    def get_secs(self):
        ''' Get all neuron sections '''
        secs = [sec for sec in h.nodeP]
        secs.extend([sec for sec in h.MYSAP])
        secs.extend([sec for sec in h.FLUTP])
        secs.extend([sec for sec in h.STINP])

        secs.extend([sec for sec in h.nodeC])
        secs.extend([sec for sec in h.MYSAC])
        secs.extend([sec for sec in h.FLUTC])
        secs.extend([sec for sec in h.STINC])

        secs.extend([sec for sec in h.nodeT])
        secs.extend([sec for sec in h.MYSAT])
        secs.extend([sec for sec in h.FLUTT])
        secs.extend([sec for sec in h.STINT])
        
        secs.extend([sec for sec in h.iseg])
        secs.extend([sec for sec in h.soma])

        return secs

    def get_nodesP(self):
        ''' Get all peripheral node sections '''
        secs = [sec for sec in h.nodeP]
        return secs

    def get_nodesC(self):
        ''' Get all central node sections '''
        secs = [sec for sec in h.nodeC]
        return secs

    def get_nodesT(self):
        ''' Get all stem node sections '''
        secs = [sec for sec in h.nodeT]
        return secs
    
    def get_mysaP(self):
        ''' Get all peripheral mysa sections '''
        return [sec for sec in h.MYSAP]

    def get_mysaC(self):
        ''' Get all central mysa sections '''
        return [sec for sec in h.MYSAC]

    def get_mysaT(self):
        ''' Get all stem mysa sections '''
        return [sec for sec in h.MYSAT]
    
    def get_flutP(self):
        ''' Get all peripheral flut sections '''
        return [sec for sec in h.FLUTP]

    def get_flutC(self):
        ''' Get all central flut sections '''
        return [sec for sec in h.FLUTC]

    def get_flutT(self):
        ''' Get all stem flut sections '''
        return [sec for sec in h.FLUTT]
    
    def get_stinP(self):
        ''' Get all peripheral stin sections '''
        return [sec for sec in h.STINP]

    def get_stinC(self):
        ''' Get all central stin sections '''
        return [sec for sec in h.STINC]

    def get_stinT(self):
        ''' Get all stem stin sections '''
        return [sec for sec in h.STINT]
    
    def get_iseg(self):
        ''' Get all iseg sections '''
        return [sec for sec in h.iseg]
    
    def get_soma(self):
        ''' Get all soma sections '''
        return [sec for sec in h.soma]

    
    def _construct_cell(self):
        """ mammalian motor axon """
        
        import numpy as np

        # set up variables

        self.axon_trajectory_iseg = self.get_variable('axon_trajectory_iseg')
        self.axon_trajectory_stem = self.get_variable('axon_trajectory_stem')
        self.axon_trajectory_per = self.get_variable('axon_trajectory_per')
        self.axon_trajectory_cen = self.get_variable('axon_trajectory_cen')
        
        self.fiberD_iseg = self.get_variable('fiberD_iseg')
        self.fiberD_stem = self.get_variable('fiberD_stem')
        self.fiberD_per = self.get_variable('fiberD_per')
        self.fiberD_cen = self.get_variable('fiberD_cen')
        
        
        # Initial Segment
        
        # set up parameters for iseg
        self.axon_trajectory_iseg = find_devor_compart_coords_high_res(self.fiberD_iseg,self.axon_trajectory_iseg,'iseg') # create 1 um spacing
        
        self.numIseg = shape(self.axon_trajectory_iseg)[0]-1 # shape-1 points to provide endpoint for last iseg compartment
        
        # find exact distance between each iseg node point
        self.iseg_indv_len = np.zeros(shape(self.axon_trajectory_iseg)[0]-1)
        for p in range(shape(self.axon_trajectory_iseg)[0]-1):
            self.iseg_indv_len[p] = sqrt(sum((self.axon_trajectory_iseg[p+1,:]-self.axon_trajectory_iseg[p,:])**2))
        
        # load iseg parameters into hoc
        h('numIseg = %i' % self.numIseg) # total number of inital segment nodes
        h('fiberD_iseg = %i' % self.fiberD_iseg) # diameter of initial segment
        h('objref iseg_indv_len') # length of individual initial segments
        h('iseg_indv_len = new Vector(%i)' % shape(self.iseg_indv_len)[0])
        for i in range(shape(self.iseg_indv_len)[0]):
            h.iseg_indv_len.x[i] = self.iseg_indv_len[i]
        
        # load iseg coordinates into hoc
        h('objref nxI')
        h('objref nyI')
        h('objref nzI')
        h('nxI = new Vector(%i)' % shape(self.axon_trajectory_iseg)[0])
        h('nyI = new Vector(%i)' % shape(self.axon_trajectory_iseg)[0])
        h('nzI = new Vector(%i)' % shape(self.axon_trajectory_iseg)[0])
        for i in range(shape(self.axon_trajectory_iseg)[0]):
            h.nxI.x[i] = self.axon_trajectory_iseg[i,0]
            h.nyI.x[i] = self.axon_trajectory_iseg[i,1]
            h.nzI.x[i] = self.axon_trajectory_iseg[i,2]

        
        # Stem Axon
        
        # Find node coordinates for lower stem axon
        self.NODE_COORDINATES_stem = find_devor_compart_coords_high_res(self.fiberD_stem,self.axon_trajectory_stem,branch='stem_with_iseg')
        
        # determine length of lower stem axon
        length = 0
        for p in range(shape(self.axon_trajectory_stem)[0]-1):
            length = length + sqrt(sum((self.axon_trajectory_stem[p+1,:]-self.axon_trajectory_stem[p,:])**2))*1e3 # convert mm to um
        self.stem_len = length
        
        # determine internode lengths
        base_inter_lengths = [86,131,169,202] # Amir & Devor 2003, added 1,1,1,1 (4 total) to account for node lengths
        total_inter_length = sum(base_inter_lengths)
        percentages = list()
        for n in base_inter_lengths:
            percentages.append(n/total_inter_length)
        self.varLenT = [n*(self.stem_len-2) for n in percentages] # leave space for final node


        # topological parameters, based on fiberD and NODE_COORDINATES (i.e. axon trajectory)
        self.axonnodes_stem = 1+shape(self.NODE_COORDINATES_stem)[0]//11
        self.paranodes1_stem = 2*(self.axonnodes_stem-1)
        self.paranodes2_stem = 2*(self.axonnodes_stem-1)
        self.stins_stem = 6*(self.axonnodes_stem-1)
        self.axontotal_stem = self.axonnodes_stem+self.paranodes1_stem+self.paranodes2_stem+self.stins_stem

        # define parameters
        h('fiberD_stem = %.1f' % self.fiberD_stem)
        h('axonnodesT = %i' % self.axonnodes_stem)
        h('paranodes1T = %i' % self.paranodes1_stem)
        h('paranodes2T = %i' % self.paranodes2_stem)
        h('stinsT = %i' % self.stins_stem)
        h('axontotalT = %i' % self.axontotal_stem)

        h('objref varLenT')
        h('varLenT = new Vector(%i)' % shape(self.varLenT)[0])
        for i in range(shape(self.varLenT)[0]):
            h.varLenT.x[i] = self.varLenT[i]

        # load node coordinates into hoc
        h('objref nxT')
        h('objref nyT')
        h('objref nzT')
        h('nxT = new Vector(%i)' % shape(self.NODE_COORDINATES_stem)[0])
        h('nyT = new Vector(%i)' % shape(self.NODE_COORDINATES_stem)[0])
        h('nzT = new Vector(%i)' % shape(self.NODE_COORDINATES_stem)[0])
        for i in range(shape(self.NODE_COORDINATES_stem)[0]):
            h.nxT.x[i] = self.NODE_COORDINATES_stem[i,0]
            h.nyT.x[i] = self.NODE_COORDINATES_stem[i,1]
            h.nzT.x[i] = self.NODE_COORDINATES_stem[i,2]
            
        
        # Peripheral Axon
        
        # Find node coordinates for peripheral axon
        self.NODE_COORDINATES_per = find_devor_compart_coords_high_res(self.fiberD_per,self.axon_trajectory_per,branch='peripheral')

        # topological parameters, based on fiberD and compartment coordinates (i.e. axon trajectory)
        self.axonnodes_per = 1+shape(self.NODE_COORDINATES_per)[0]//11 # truncated to account for only node coordinates
        # determine number of mysa, cut off at last node
        self.paranodes1_per = 2*(self.axonnodes_per-1)
        # determine number of flut
        self.paranodes2_per = 2*(self.axonnodes_per-1)
        # determine number of stin
        self.stins_per = 6*(self.axonnodes_per-1)
        # get total number of compartments
        self.axontotal_per = self.axonnodes_per+self.paranodes1_per+self.paranodes2_per+self.stins_per

        # define parameters
        h('fiberD_peripheral = %.1f' % self.fiberD_per)
        h('axonnodesP = %i' % self.axonnodes_per)
        h('paranodes1P = %i' % self.paranodes1_per)
        h('paranodes2P = %i' % self.paranodes2_per)
        h('stinsP = %i' % self.stins_per)
        h('axontotalP = %i' % self.axontotal_per)

        # load node coordinates into hoc
        h('objref nxP')
        h('objref nyP')
        h('objref nzP')
        h('nxP = new Vector(%i)' % shape(self.NODE_COORDINATES_per)[0])
        h('nyP = new Vector(%i)' % shape(self.NODE_COORDINATES_per)[0])
        h('nzP = new Vector(%i)' % shape(self.NODE_COORDINATES_per)[0])
        for i in range(shape(self.NODE_COORDINATES_per)[0]):
            h.nxP.x[i] = self.NODE_COORDINATES_per[i,0]
            h.nyP.x[i] = self.NODE_COORDINATES_per[i,1]
            h.nzP.x[i] = self.NODE_COORDINATES_per[i,2]
            
            
        # Central Axon
        
        # Find node coordinates for central axon
        self.NODE_COORDINATES_cen = find_devor_compart_coords_high_res(self.fiberD_cen,self.axon_trajectory_cen,branch='central')

        # topological parameters, based on fiberD and NODE_COORDINATES (i.e. axon trajectory)
        self.axonnodes_cen = 1+shape(self.NODE_COORDINATES_cen)[0]//11 # truncated to account for only node coordinates
        # determine number of mysa, cut off at last node
        self.paranodes1_cen = 2*(self.axonnodes_cen-1)
        # determine number of flut
        self.paranodes2_cen = 2*(self.axonnodes_cen-1)
        # determine number of stin
        self.stins_cen = 6*(self.axonnodes_cen-1)
        # get total number of compartments
        self.axontotal_cen = self.axonnodes_cen+self.paranodes1_cen+self.paranodes2_cen+self.stins_cen

        # define parameters
        h('fiberD_central = %.1f' % self.fiberD_cen)
        h('axonnodesC = %i' % self.axonnodes_cen)
        h('paranodes1C = %i' % self.paranodes1_cen)
        h('paranodes2C = %i' % self.paranodes2_cen)
        h('stinsC = %i' % self.stins_cen)
        h('axontotalC = %i' % self.axontotal_cen)

        # load node coordinates into hoc
        h('objref nxC')
        h('objref nyC')
        h('objref nzC')
        h('nxC = new Vector(%i)' % shape(self.NODE_COORDINATES_cen)[0])
        h('nyC = new Vector(%i)' % shape(self.NODE_COORDINATES_cen)[0])
        h('nzC = new Vector(%i)' % shape(self.NODE_COORDINATES_cen)[0])
        for i in range(shape(self.NODE_COORDINATES_cen)[0]):
            h.nxC.x[i] = self.NODE_COORDINATES_cen[i,0]
            h.nyC.x[i] = self.NODE_COORDINATES_cen[i,1]
            h.nzC.x[i] = self.NODE_COORDINATES_cen[i,2]
        
        
        # Multicompartment Soma
        somaL = 40
        somaD = 40
        h('somaL = %i' % somaL)
        h('somaD = %i' % somaD)
        loc = np.arange(0.5,somaL,1) # set up locations along soma length
        rad = []
        for x in loc:
            rad.append(np.sqrt((somaD/2)**2 - (x - (somaL/2))**2)) # calculate radius for each soma segment location
        h('objref somadiams') # load somadiams vector into hoc
        h('somadiams = new Vector(%i)' % shape(rad)[0])
        for i in range(shape(rad)[0]):
            h.somadiams.x[i] = 2*rad[i]

        # Load hoc file
        h.load_file(self.CELL_DIR+self.CELL_NAME)
        self.root = h.nodeT[0]
--------------------------------------------------------------------------------------------------
DRG A-DELTA MICROELECTRODE HYBRID RECORDING MODEL
Lauren R. Madden, Robert D. Graham, Scott F. Lempka, Tim M. Bruns | University of Michigan | 2023
README.txt

Manuscript: L.R. Madden, R.D. Graham, S.F. Lempka, T.M. Bruns, "Multiformity of
extracellular microelectrode recordings from Aδ neurons in the dorsal root ganglia: A
computational modeling study." 2023.

We built a computational model of A-delta LTMR neurons in the dorsal root ganglia to investigate
factors that influence the extracellular waveforms recorded by microelectrodes. Our model
demonstrates how the unique structure of these neurons can lead to diverse and often multiphasic
waveform profiles depending on the location of the recording contact relative to microanatomical
neural components. Our model also provides insight into the neurophysiological function of axon
glomeruli that are often present in these neurons.

This code allows the user to generate simulated extracellular (and intracellular) voltage
waveforms of A-delta LTMR NEURON models with varying glomerular initial segment (GIS)
trajectory types, GIS lengths, and microelectrode array positions. Requires installation of
NEURON simulation environment (https://neuron.yale.edu/neuron/). Simulation of additional GIS
trajectories and microelectrode array positions also requires licenses for both COMSOL
Multiphysics® software and MATLAB by MathWorks.

Note: the directory NEURON/ADelta_files with the .hoc file and .mod files must be compiled via
mknrndll after downloading the files. The resulting .dll file will be put into the same
ADelta_files directory, and must be kept there.
Additionally, this repository's organizational structure needs to be maintained in order for
simulations to run properly.
--------------------------------------------------------------------------------------------------

Created with NEURON Version 8.2.0 and Python Version 3.7


------ Available via NEURON with Python installation only ------

--- Run Simulation with pre-generated interpolations
	- Provided pre-generated interpolations:
		- Grey matter/extracellular conductivity (sigma_ex) = 0.26 S/m
		- Ramon y Cajal 1911 tracing of glomerular initial segment (GIS) trajectory
		- Found in InterpolationFiles directory:
			- interpolation_e36_sigma0p26.txt
			- interpolation_e38_sigma0p26.txt
			- interpolation_e40_sigma0p26.txt
			- interpolation_e42_sigma0p26.txt
			- interpolation_e44_sigma0p26.txt
			- interpolation_e46_sigma0p26.txt
	- Run bash shell.sh for simulations with standard protocol
		- Generates plot found in 15 um panels in Figure 9 A and C
	- Run bash shell_withAEC.sh for simulations with auto-ephaptic coupling (AEC)
		- Note: each simulation with auto-ephaptic coupling can take 2-3 hours



------ Requires both MATLAB and COMSOL licenses ------
Created with MATLAB Version R2020b and COMSOL Version 6.0

--- 1. COMSOL FEM Model Solutions
	The directory FEM/COMSOL_files contains base model files that can be meshed and solved
with COMSOL-MATLAB Simulink scripts found in FEM/MATLAB-COMSOL-Simulink_scripts. Requires both
MATLAB and COMSOL licenses.
	- arrayLocationShift.m shifts the location of the microelectrode array and changes the
		active electrode, and generates model files with solutions
			- Corresponds with Figure 9 
	- electrode_and_conductivity_change.m solves for different grey matter conductivities
		and active electrodes, and generates model files with solutions
			- Corresponds with Figures 5 and 6

--- 2. Generate NEURON Model Coordinate Files
	- NEURON/XYZcoords contains coordinate files for each axon branch
	- To generate a new GIS trajectory (i.e. iSeg_X.txt), use one of the following:
		- glomerulus_geom_spiral.py
			- Generates a spiral trajectory given tortuosity (and other parameters)
			- Corresponds with Figure 5
		- glomerlus_geom_spinral_same_height_dif_len.py
			- Generates a spiral trajectory given length (and other parameters)
			- Corresponds with Figure 8
		- glomerulus_geom_randwalk.py
			- Generates a random walk trajectory given tortuosity (and other
				parameters)
			- Corresponds with Figure 6
			- Note: this algorithm generates multiple iterations of a path, and
				can sometimes get stuck. Just quit and re-run the script
				when this happens.
		- Note: edit these .py files to adjust parameters
	- Next, generate the XYZ coordinates for the new NEURON model using either:
		- generate_nrn_centers_tracing.py
			- generates a nrn_centers.txt file using the original tracing coordinates
				(as found in NEURON/XYZcoords/iSeg.txt)
		- generate_nrn_centers_variedGIS.py
			- generates a nrn_centers_X.txt file using a previously generated and
				specified glomerulus trajectory file (iSeg_0_X.txt) and creates
				a new file adjusted for global positioning (iSeg_X.txt)
			- Note: edit this .py file to specify GIS trajectory suffix

--- 3. Interpolate NEURON Model with FEM Model
	- Use MATLAB Simulink Scripts in FEM/MATLAB-COMSOL-Simulink_scripts
		- interpolate_and_export_arrayShift.m interpolates for different microelectrode
			positions and active contacts. Edit the nrn_centers.txt file name for a
			newly generated GIS trajectory (nrn_centers_X.txt).
		- interpolate_and_export_e_and_sigma.m interpolates for different grey matter
			conductivity values and contacts. Edit the nrn_centers.txt file name
			for a newly generated GIS trajectory (nrn_centers_X.txt).

--- 4. Run NEURON Simulation
	- Edit the shell_params.sh script to indicate the iSeg_X.txt file and corresponding
		interpolation_X.txt file
	- To simulate with the standard protocol, run bash shell_params.sh
	- To simulate with auto-ephaptic coupling, run bash shell_params_withAEC.sh


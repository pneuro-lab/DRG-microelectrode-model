#!/bin/bash

# Run Simulation with given interpolation.txt file and iSeg.txt file, standard protocol
# Example using pre-made iSeg trajectory and FEM interpolation is given

# exracellular voltage recording
python runSimulation.py "extra" "iSeg.txt" "interpolation_e36_sigma0p26.txt"
python view_results_params.py

# intracellular voltage recording (does not require interpolation file)
python runSimulation.py "intra" "iSeg.txt"
python view_results_params.py

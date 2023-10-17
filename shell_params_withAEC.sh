#!/bin/bash

# Run Simulation with given interpolation.txt file and iSeg.txt file, uses auto-ephaptic coupling
# Example using pre-made iSeg trajectory and FEM interpolation is given

# extracellular voltage recording
python runSimulation_withAEC.py "extra" "iSeg.txt" "interpolation_e36_sigma0p26.txt"
python view_results.py

# intracellular voltage recording (does not require interpolation file)
# requires extracellular conductivity input in S/m
python runSimulation_withAEC.py "intra" "iSeg.txt" 0.26
python view_results.py
#!/bin/bash

# Run Simulation with pre-made interpolations, uses auto-ephaptic coupling

# extracellular voltage recordings
python runSimulation_withAEC.py "extra" "iSeg.txt" "interpolation_e36_sigma0p26.txt"
python runSimulation_withAEC.py "extra" "iSeg.txt" "interpolation_e38_sigma0p26.txt"
python runSimulation_withAEC.py "extra" "iSeg.txt" "interpolation_e40_sigma0p26.txt"
python runSimulation_withAEC.py "extra" "iSeg.txt" "interpolation_e42_sigma0p26.txt"
python runSimulation_withAEC.py "extra" "iSeg.txt" "interpolation_e44_sigma0p26.txt"
python runSimulation_withAEC.py "extra" "iSeg.txt" "interpolation_e46_sigma0p26.txt"
python view_results.py

# intracellular voltage recording (does not require interpolation file)
# requires extracellular conductivity input in S/m
python runSimulation_withAEC.py "intra" "iSeg.txt" 0.26
python view_results.py
#!/bin/bash

# Run Simulation with pre-made interpolations, standard protocol

# extracellular voltage recordings
python runSimulation.py "extra" "iSeg.txt" "interpolation_e36_sigma0p26.txt"
python runSimulation.py "extra" "iSeg.txt" "interpolation_e38_sigma0p26.txt"
python runSimulation.py "extra" "iSeg.txt" "interpolation_e40_sigma0p26.txt"
python runSimulation.py "extra" "iSeg.txt" "interpolation_e42_sigma0p26.txt"
python runSimulation.py "extra" "iSeg.txt" "interpolation_e44_sigma0p26.txt"
python runSimulation.py "extra" "iSeg.txt" "interpolation_e46_sigma0p26.txt"
python view_results.py

# intracellular voltage recording (does not require interpolation file)
python runSimulation.py "intra" "iSeg.txt"
python view_results.py
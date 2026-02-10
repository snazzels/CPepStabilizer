#!/usr/bin/env python3
"""
## System Setup
"""
import numpy as np
import mdtraj as md
import pandas as pd
from functions import *

# Open the results.txt file to capture output from print statements
results_file = open('results.txt', 'w')

"""
## Calculate the total Free Binding Energies per Compound and Run
"""
AL_gbsa_path = 'AL_out/AL_output.dat'
BL_gbsa_path = 'BL_out/BL_output.dat'
ABL_gbsa_path = 'AB_L_out/AB_L_output.dat'

AL_dG, AL_sdG, AL_seG = read_output(AL_gbsa_path, gbsa_pattern)
BL_dG, BL_sdG, BL_seG = read_output(BL_gbsa_path, gbsa_pattern)
ABL_dG, ABL_sdG, ABL_seG = read_output(ABL_gbsa_path, gbsa_pattern)

print(f'dG_ABL = {ABL_dG} +- {ABL_sdG}', file=results_file)
print(f'dGAL = {AL_dG} +- {AL_sdG}', file=results_file)
print(f'dGBL = {BL_dG} +- {BL_sdG}', file=results_file)


# Close the file after all print statements
results_file.close()


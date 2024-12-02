#!/usr/bin/env python3
import os
import sys
import re
import glob
import math
fastcar_dir = "/home/hgeindre/ProjectFASTCAR/FASTCAR2"
sys.path.append(fastcar_dir)
import readparam as rp
import readoutput as ro
import filemanager as fm
import misc

# This script extracts data from final calculations 
# It will calculate local electrophilicity if NPA selected
# It will print BD* in case of NBO
# In other cases it will extract energies

# Get information from parameters_mod.txt file
paramfile = rp.Parameters('../parameters.tmp')
sumfilename = "summary.log"

Energies = []
NPA = []
NBO = []
local_E = []

for filename in sorted(glob.glob(f'*{paramfile.experience_number}*.out')):
    outputfile = ro.Output(filename)
    outputfile.readnrj()
    # Extract NPA
    if 'cdft' in filename:
        outputfile.readcharges()
    # Extract NBO
    elif 'nbo' in filename:
        outputfile.readnbo()
        NBO = outputfile.NBO
    # Extract energies from other calculations
    else:   
        outputfile.readallnrjs()
        Energies.append(outputfile.allnrjs)

# Write in a new file 'results' all the data extracted
fm.writeresults(Energies,NPA,local_E,NBO)
with open(f'../{sumfilename}', 'a') as sumfile:
    sumfile.write('Informations about final calculations written in '
                  'results.txt\n')
    ending = 'Destination reached'
    sumfile.write(f'{ending:-^80}\n')
#os.remove('../parameters.tmp')

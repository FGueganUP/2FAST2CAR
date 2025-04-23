#!/usr/bin/env python3
import os
import sys
import datetime
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

sumfilename = "summary.log"

now = datetime.datetime.now()
date = now.strftime("%d/%m/%y %H:%M")
with open(f'../{sumfilename}','a') as sumfile:
    sumfile.write(f'Calling wheels.py {date}\n')

# Get information from parameters_mod.txt file
paramfile = rp.Parameters('../parameters.tmp')

energies = []
NPA = []
NBO = []
local_E = []

#for filename in sorted(glob.glob(f'*{paramfile.calc_name}*.out')):
#    outputfile = ro.Output(filename)
#    # Extract NBO
#    if 'nbo' in filename:
#        outputfile.readnbo()
#        NBO = outputfile.NBO
#    # Extract energies from other calculations
#    else:   
#        energies.append(outputfile.energies)

# Write in a new file 'results' all the data extracted
#fm.writeresults(energies,NPA,local_E,NBO)

nrjs = fm.readsummary('../structures.log')

with open(f'../{sumfilename}', 'a') as sumfile:
    sumfile.write('Informations about final calculations written in '
                  'results.txt\n')
    sumfile.write(f'{len(nrjs)} unique conformers found and optimised\n')
    ending = 'Destination reached'
    sumfile.write(f'{ending:-^80}\n')
#os.remove('../parameters.tmp')

#!/usr/bin/env python3
import sys
import os
import glob
import re
import shutil
fastcar_dir = "/home/hgeindre/ProjectFASTCAR/FASTCAR2"
sys.path.append(fastcar_dir)
import readoutput as ro
import readparam as rp
import filemanager as fm
#
# This script is intended to take as single argument a calculation in Gaussian 16:
#
#  - for TS a frequency calculation is the minimum criteria 
#  - for geometry an optimisation that converged is the minimum criteria
#
# The parameters for the calculations must be given in a parameters.txt file

sumfilename = "summary.log"

# Get the file name from the command-line argument
if len(sys.argv) == 1:
    if os.path.isdir('Final/'):
        file_name = glob.glob('Final/*_lowest.out')[0]
        wora = 'a'
        loop = True
    else:
        for line in piclines:
            print(line.rstrip())
        sys.exit(1)
elif len(sys.argv) != 2:
    print("Please give one and only one argument, the name of an ouputfile")
    sys.exit(1)
else:
    file_name = sys.argv[1]
    loop = False
    wora = 'w'

# Open output file to extract informations
startoutput = ro.Output(file_name)

# Read frequencies and check if the first one is imaginary
startoutput.readfreqs() 
# Read energy
startoutput.readnrj() 

# Get information from parameters.txt file
if loop:
    paramfile = rp.Parameters('parameters.tmp',startoutput)
    paramfile.check_st += 1
else:
    paramfile = rp.Parameters('parameters.txt',startoutput)
    paramfile.check_st = 0

with open(sumfilename,wora) as sumfile:
    emp = ''
    if wora == 'w':
        heading = 'FASTCAR calculation summary'
        sumfile.write(f'{emp:-^80}\n')
        sumfile.write(f'{heading:-^80}\n')
        sumfile.write(f'{emp:-^80}\n\n\n')
    heading = f'Loop {paramfile.check_st}'
    sumfile.write(f'{heading:-^80}\n\n')
    heading = 'Starting point'
    sumfile.write(f'{heading:-^80}\n')
    sumfile.write(f'Starting file: {startoutput.file_name}\n')
    sumfile.write(f'Software: {startoutput.software}\n')
    sumfile.write(f'Normal termination: {startoutput.normterm}\n')
    sumfile.write(f'Charge: {startoutput.charge}, multiplicity: '
                  f'{startoutput.multiplicity}\n')
    sumfile.write(f'Number of atoms: {startoutput.num_atom}\n')
    if startoutput.ts:
        struc = "transition state"
    else:
        struc = "minimum"
    sumfile.write(f'Structure type: {struc}\n')
    sumfile.write(f'Energy: {startoutput.last_SCF} ua\n')
    if startoutput.ts:
        sumfile.write(f'Im. Freq.: {startoutput.freqs[0]} cm-1\n')
    sumfile.write(f'{emp:-^80}\n\n')
    heading = 'CREST parameters'
    sumfile.write(f'{heading:-^80}\n')
    sumfile.write(f'Crest version: {paramfile.version_crest}\n')
    sumfile.write(f'Solvent: {paramfile.solvent_crest}\n')
    sumfile.write(f'E win: {paramfile.ewin_crest}kcal/mol\n')
    sumfile.write(f'NCI: {paramfile.nci}\n')
    if startoutput.ts:
        sumfile.write(f'constraints:\n')
        for bonds in paramfile.bonds_to_freeze:
            sumfile.write(f'  distance: {bonds[0]}, {bonds[1]}, {bonds[2]} \n')
        for angles in paramfile.angles_to_freeze:
            sumfile.write(f'  angle: {angles[0]}, {angles[1]}, {angles[2]}, '
                          f'{angles[3]} \n')
        for dihedral_angles in paramfile.dihedrals_to_freeze:
            sumfile.write(f'  dihedral: {dihedral_angles[0]}, '
                          f'{dihedral_angles[1]}, {dihedral_angles[2]}, '
                          f'{dihedral_angles[3]}, {dihedral_angles[4]} \n')
    sumfile.write(f'{emp:-^80}\n\n')


# Create struc.xyz for CREST calculation
with open('struc.xyz', 'w') as file:
    ecs = fm.readepc()
    file.write(f'{startoutput.num_atom}\n')
    file.write('\n')
    for c,coord in enumerate(startoutput.coordinates):
        file.write(f'{ecs[startoutput.atnums[c]][0]} ')
        file.write(' '.join(map(str, coord)))
        file.write('\n')

# Construction of constraints.inp for CREST calculation
if startoutput.ts:
    fm.writeconstraints(paramfile)

# File with all parameters chosen
fm.writeparam(startoutput,paramfile)

# Construction of script.sub to launch CREST calculation and call engine.py
fm.writecrestsubscript(startoutput,paramfile)

# Construction of CREST folder
if paramfile.check_st:
    crest_dir = f'CREST_{paramfile.check_st}'
else:
    crest_dir = 'CREST'
if not paramfile.woc:
    if os.path.isdir(crest_dir):
        shutil.rmtree(crest_dir)
    os.mkdir(crest_dir)
shutil.move('struc.xyz', f'{crest_dir}/struc.xyz')
shutil.move('script_CREST.sub', f'{crest_dir}/script_CREST.sub')
if startoutput.ts:
    shutil.move('constraints.inp',f'{crest_dir}/constraints.inp')

# Launch calculation
os.chdir(crest_dir)
os.system(f'sbatch script_CREST.sub')
with open(f'../{sumfilename}','a') as sumfile:
    sumfile.write(f'CREST calculation submitted in {crest_dir}/ directory, '
                   'results in CrestAnalysis.txt\n') 


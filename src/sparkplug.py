#!/usr/bin/env python3
import sys
import os
import argparse
import glob
import re
import shutil
import datetime
fastcar_dir = "/home/hgeindre/ProjectFASTCAR/FASTCAR2"
sys.path.append(fastcar_dir)
import readoutput as ro
import readparam as rp
import readconfig as rc
import filemanager as fm
#
# This script is used to start a new fastcar calculation, or start a new loop
# The parameters for the calculations must be given in a parameters.txt file

# Some options that can be specified when launching a calculation
parser = argparse.ArgumentParser()
parser.add_argument("filename", type=str, nargs='?',
                    help="The name of the output file to start the calculation"
                    "no name will start a new loop or show a surprise")
parser.add_argument("-ac", "--addcalc", action="store_true", 
                    default=False, help="Just perform additional calculations")
args = parser.parse_args()

sumfilename = "summary.log"

# Only perform additional calculation
if args.addcalc:
    if not os.path.isdir('Final/'):
        print('No Final/ directory found, a complete calculation must already '
              'have been done to perform only the additional calculations')
        sys.exit(1)
    with open(f'{sumfilename}','a') as sumfile:
        emp = ''
        heading = 'New additional calculations'
        sumfile.write(f'\n\n{heading:-^80}\n')
    file_name = glob.glob('Final/*_lowest.out')[0]
    startoutput = ro.Output(file_name)
    paramfile = rp.Parameters('parameters.txt')
    tmp_paramfile = rp.Parameters('parameters.tmp')
    clust = rc.Cluster(fastcar_dir)
    fm.updateparam(startoutput,paramfile,tmp_paramfile,clust)
    paramfile = rp.Parameters('parameters.tmp')
    os.chdir('Final/')
    scriptname = fm.writesubscript('abs',paramfile,clust)
    os.system(f"{clust.subco} {scriptname}")   
    sys.exit(0)
     
# Get the file name from the command-line argument, if no argument either loop
# if Final directory already exists or show a nice ASCII art
if len(sys.argv) == 1:
    if os.path.isdir('Final/'):
        file_name = glob.glob('Final/*_lowest.out')[0]
        wora = 'a' #Write or Append in the summary file
        loop = True
    else:
        with open(f'{fastcar_dir}/pic.txt','r') as picfile:
            piclines = picfile.readlines()
            for line in piclines:
                print(line.rstrip())
            sys.exit(0)
elif len(sys.argv) != 2:
    print("Please give one and only one argument, the name of an ouputfile")
    sys.exit(1)
else:
    file_name = args.filename
    loop = False
    wora = 'w' #Write or Append in the summary file
    # Delete all previous CREST directories if it's a new calculation
    crest_dirs = [f.name for f in os.scandir('.') if f.is_dir() \
                  and f.name.startswith('CREST')]
    for cd in crest_dirs:
        shutil.rmtree(cd)

# Write some informations in the summary file
now = datetime.datetime.now()
date = now.strftime("%d/%m/%y %H:%M")
with open(f'{sumfilename}',wora) as sumfile:
    if wora == 'w':
        emp = ''
        heading = 'FASTCAR calculation summary'
        sumfile.write(f'{emp:-^80}\n')
        sumfile.write(f'{heading:-^80}\n')
        sumfile.write(f'{emp:-^80}\n\n')
    sumfile.write(f'Calling sparkplug.py {date}\n\n')

# Open output file to extract informations
startoutput = ro.Output(file_name)

# Read informations about cluster in file config.txt
clust = rc.Cluster(fastcar_dir)

# The file provided must have a freq calculation
if not hasattr(startoutput,'freqs'):
    print("Please provide an output file with a frequency calculation")
    sys.exit(1)

# Get information from parameters.txt or parameters.tmp file
if loop:
    paramfile = rp.Parameters('parameters.tmp')
    paramfile.getgeneral()
    paramfile.getwoc()
    paramfile.loopn += 1
else:
    paramfile = rp.Parameters('parameters.txt',startoutput,clust)
    paramfile.getopt()
    paramfile.getwoc()
    if paramfile.software == 'same':
        for soft_options in clust.default:
            if soft_options in startoutput.software.lower():
                paramfile.software = clust.default[soft_options]
                break
        else:
            print(f"There is no default version for the software "
                   "{startoutput.software} in the file config.txt")
            sys.exit(1)
    paramfile.loopn = 0
    #Create directory for opt calculations  
    if startoutput.ts:
        calcdir = 'TS-CREST'
        calctype = 'opt-ts'
    else:
        calcdir = 'Geo-CREST'
        calctype = 'opt'

    if os.path.isdir(f'{calcdir}'):
        if not paramfile.woc:
            shutil.rmtree(f'{calcdir}')
            os.mkdir(f'{calcdir}')
    else:
        os.mkdir(f'{calcdir}')

    # Check if given outputfile as the same level of theory than the one that
    # will be used for opt calculations, if yes the outputfile is included
    # in the results, if not it is reoptimised
    jobname = f"{startoutput.base_name}-0"
    reopt = fm.comptheolvl(startoutput,paramfile)
    if reopt:
        kwds = rc.Keywords(fastcar_dir)
        os.chdir(f'{calcdir}')
        inpname = f'{jobname}.inp' 
        paramfile.charge = startoutput.charge
        paramfile.mult = startoutput.mult
        kw = (kwds.keywords[paramfile.software][calctype],calctype)
        fm.writeinp(inpname,kw,startoutput.coordinates,startoutput.atnums,
                    paramfile)
        fm.writecalcsubscript(paramfile,jobname,clust)
        if not paramfile.woc:
            os.system(f"{clust.subco} {jobname}.sub")
        os.chdir('../')
        with open(f'{sumfilename}','a') as sumfile:
            sentence = (f'Starting file {startoutput.file_name} '
                       f'({startoutput.method}/{startoutput.basis_set} '
                       f'solvent = {startoutput.solvent}) reoptimised at the'
                        ' same level of theory specified in the parameters '
                       f'file ({paramfile.method}/{paramfile.basis_set} '
                       f'{paramfile.solvent}) in {calcdir} as {inpname}.\n')
            #Write only on 80 columns max
            spaces = [p for p,c in enumerate(sentence) if c == ' ']
            start = 0
            for n,sp in enumerate(spaces):
                if start+80 > len(sentence):
                    sumfile.write(f'{sentence[start:].strip()}\n\n')
                    break
                if sp > start+80:
                    sumfile.write(f'{sentence[start:spaces[n-1]].strip()}\n')
                    start = spaces[n-1]
    else:
        outname = f'{jobname}.out' 
        shutil.copy(file_name,f'{calcdir}/{outname}') 
        with open(f'{sumfilename}','a') as sumfile:
            sentence = (f'Starting file {startoutput.file_name} '
                       f'({startoutput.method}/{startoutput.basis_set} '
                       f'solvent = {startoutput.solvent}) already at the '
                        'same level of theory specified in the parameters '
                       f'file ({paramfile.method}/{paramfile.basis_set} '
                       f'solvent = {paramfile.solvent}), copied in {calcdir}'
                       f' as {outname}.\n')
            spaces = [p for p,c in enumerate(sentence) if c == ' ']
            start = 0
            #Write only on 80 columns max
            for n,sp in enumerate(spaces):
                if start+80 > len(sentence):
                    sumfile.write(f'{sentence[start:].strip()}\n\n')
                    break
                if sp > start+80:
                    sumfile.write(f'{sentence[start:spaces[n-1]].strip()}\n')
                    start = spaces[n-1]

# Get CREST parameters and constraints if TS from parameters file
paramfile.getcrest()
if startoutput.ts:
    paramfile.getconstraints(startoutput)

# Write information about starting point and parameters in summary file
with open(sumfilename,'a') as sumfile:
    emp = ''
    heading = f'Loop {paramfile.loopn}'
    sumfile.write(f'{heading:~^80}\n\n')
    heading = 'CREST calculation'
    sumfile.write(f'{emp:-^80}\n')
    sumfile.write(f'{heading:-^80}\n\n')
    heading = 'Starting point'
    sumfile.write(f'{heading:-^80}\n')
    sumfile.write(f'Starting file: {startoutput.file_name}\n')
    sumfile.write(f'Software: {startoutput.software}\n')
    sumfile.write(f'Charge: {startoutput.charge}, multiplicity: '
                  f'{startoutput.mult}\n')
    sumfile.write(f'Number of atoms: {startoutput.num_atom}\n')
    if startoutput.ts:
        struc = "transition state"
    else:
        struc = "minimum"
    sumfile.write(f'Structure type: {struc}\n')
    sumfile.write(f'Energy: {startoutput.last_SCF} ua\n')
    if startoutput.ts:
        sumfile.write(f'Im. Freq.: {startoutput.freqs[0]} cm-1\n\n')
    heading = 'CREST parameters'
    sumfile.write(f'{heading:-^80}\n')
    sumfile.write(f'Crest version: {paramfile.crest_version}\n')
    if paramfile.solvent_crest != 'none':
        sumfile.write(f'Solvent: {paramfile.solvent_crest}\n')
    sumfile.write(f'E win: {paramfile.ewin_crest} kcal/mol\n')
    if paramfile.nci != '':
        sumfile.write(f'NCI: {paramfile.nci}\n')
    if startoutput.ts:
        sumfile.write(f'constraints:\n')
        if paramfile.bonds:
            for bond in paramfile.bonds_to_freeze:
                sumfile.write(f'  distance: {bond[0]}, {bond[1]}, {bond[2]}\n')
        if paramfile.angles:
            for angle in parfile.angles_to_freeze:
                sumfile.write(f'  angle: {angle[0]}, {angle[1]}, {angle[2]}, '
                              f'{angle[3]} \n')
        if paramfile.dihedrals:
            for dihedral_angle in paramfile.dihedrals_to_freeze:
                sumfile.write(f'  dihedral: {dihedral_angle[0]}, '
                              f'{dihedral_angle[1]}, {dihedral_angle[2]}, '
                              f'{dihedral_angle[3]}, {dihedral_angle[4]} \n')
    sumfile.write(f'\n')


# Create struc.xyz for CREST calculation
with open('struc.xyz', 'w') as file:
    ecs = fm.readepc()
    file.write(f'{startoutput.num_atom}\n')
    file.write('\n')
    for c,coord in enumerate(startoutput.coordinates):
        if not startoutput.atnums[c] in ecs:
            print('Warning, no information on element number '
                 f'{startoutput.atnums[c]} in elementPlotCharac.data, please '
                  'add them in the file')
            sys.exit(1)
        file.write(f'{ecs[startoutput.atnums[c]][0]} ')
        file.write(' '.join(map(str, coord)))
        file.write('\n')

# Construction of constraints.inp for CREST calculation
if startoutput.ts:
    fm.writeconstraints(paramfile)

# Create or update parameters.tmp file
fm.writeparam(startoutput,paramfile,clust)

# Construction of script.sub to launch CREST calculation and call engine.py
if loop:
    crestjobname = f'crest{paramfile.loopn}_{startoutput.base_name}'
else:
    crestjobname= f'crest_{startoutput.base_name}'
fm.writecrestsubscript(paramfile,startoutput,crestjobname,clust,'engine.py')

# Construction of CREST folder
if paramfile.loopn:
    crest_dir = f'CREST_{paramfile.loopn}'
else:
    crest_dir = 'CREST'
if not paramfile.woc:
    if os.path.isdir(crest_dir):
        shutil.rmtree(crest_dir)
    os.mkdir(crest_dir)
shutil.move('struc.xyz', f'{crest_dir}/struc.xyz')
shutil.move(f'{crestjobname}.sub', f'{crest_dir}/{crestjobname}.sub')
if startoutput.ts:
    shutil.move('constraints.inp',f'{crest_dir}/constraints.inp')

# Launch CREST calculation and call engine.py
os.chdir(crest_dir)
os.system(f'{clust.subco} {crestjobname}.sub')
with open(f'../{sumfilename}','a') as sumfile:
    sumfile.write(f'CREST calculation submitted in {crest_dir}/ directory, '
                   'results in CrestAnalysis.txt\n') 


#!/usr/bin/env python3
import os
import sys
import glob
fastcar_dir = "/home/hgeindre/ProjectFASTCAR/FASTCAR2"
sys.path.append(fastcar_dir)
import readparam as rp
import readoutput as ro
import filemanager as fm
import readconfig as rc
import misc

sumfilename = "summary.log"
with open(sumfilename,'a') as sumfile:
    sumfile.write('\n')

# Get cluster specifications
clust = rc.Cluster(fastcar_dir)
# Get keywords
kwds = rc.Keywords(fastcar_dir)
# Get information from parameters.txt file
paramfile = rp.Parameters('../parameters.tmp')

lowestname = glob.glob('*_lowest.out')
file_name = lowestname[0]
lowestfile = ro.Output(file_name) 

#Create input for frequency calculation with the most stable geometry
soft = 'gaussian'
kw = kwds.keywords[soft]['freq']
base_name = file_name.split('.')[0]
fm.writeg16inp(f'freq_{base_name}.inp',paramfile,12,lowestfile.coordinates,
               lowestfile.atnums,kw)
jobname = f'freq_{base_name}'
fm.writeg16subscript(paramfile,jobname,clust)
if not paramfile.woc:
    os.system(f"sbatch {jobname}g16.sub")
with open(f'../{sumfilename}','a') as sumfile:
    sumfile.write(f'Freq calculation submitted for {file_name}\n')

# Launch IRC calculation with the most suitable TS if chosen
if paramfile.choice == 'irc':
    lowestfile.extractcoord()
    IRC = ["forward", "reverse"]
    for item in IRC:
        inpnameirc = f'{str(item)}_IRC_{paramfile.experience_number}.inp'
        soft = 'gaussian'
        kw = kwds.keywords[soft]['irc']
        kw = kw.replace('direction',item)
        kw2 = kwds.keywords[soft]['irc-opt']
        chkn = f'{paramfile.experience_number}_{str(item)}.chk'
        fm.writeg16inp(inpnameirc,paramfile,12,lowestfile.coordinates,
                       lowestfile.atnums,kw,chkname=chkn,link1=True,
                       keywords2=kw2)
        ircjobname = f'{str(item)}_IRC_{paramfile.experience_number}'
        fm.writeg16subscript(paramfile,ircjobname,clust)
        if not paramfile.woc:
            os.system(f"sbatch {ircjobname}g16.sub")
    with open(f'../{sumfilename}','a') as sumfile:
        sumfile.write(f'IRC calculations submitted for {file_name}\n')


if paramfile.choice == 'cdft':
    # Launch conceptual DFT calculations
    for i in range(1, 4): 
        inpname = f'{paramfile.choice}_{paramfile.experience_number}_{i}'+\
                  '.inp'
        chcdft = int(paramfile.charge) + i - 2
        multcdft = int(paramfile.multiplicity) + i % 2
        soft = 'gaussian'
        kw = kwds.keywords[soft]['cdft']
        fm.writeg16inp(inpname,paramfile,12,lowestfile.coordinates,
                       lowestfile.atnums,kw,chcdft,multcdft)
        jobnamecdft=f'{paramfile.choice}_{paramfile.experience_number}_{i}'
        fm.writeg16subscript(paramfile,jobnamecdft,clust)
        if not paramfile.woc:
            os.system(f"sbatch {jobnamecdft}g16.sub")
    with open(f'../{sumfilename}','a') as sumfile:
        sumfile.write(f'CDFT calculations submitted for {file_name}\n')
elif paramfile.choice == 'nbo':
    # Launch NBO calculation
    soft = 'gaussian'
    kw = kwds.keywords[soft]['nbo']
    fm.writeg16inp(inpname,paramfile,12,lowestfile.coordinates,
                   lowestfile.atnums,"pop=NBO7")
    jobnameNBO = f'nbo_{paramfile.experience_number}'
    fm.writeg16subscript(paramfile,jobnameNBO,clust,g16C=True)
    if not paramfile.woc:
        os.system(f"sbatch {jobnameNBO}g16.sub")
    with open(f'../{sumfilename}','a') as sumfile:
        sumfile.write(f'NBO calculations submitted for {file_name}\n')

# Call wheels.py
fm.writefinalscript(paramfile)
os.system(f"sbatch script_wheels.sub")


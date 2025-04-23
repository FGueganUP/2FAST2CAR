#!/usr/bin/env python3
import sys
import os
import subprocess
import shutil
import re
import glob
import datetime
import spyrmsd
from spyrmsd import io, rmsd
fastcar_dir = "/home/hgeindre/ProjectFASTCAR/FASTCAR2"
sys.path.append(fastcar_dir)
import readparam as rp
import readconfig as rc
import filemanager as fm

sumfilename = "summary.log"
uconfilename = "crest_uconfos.xyz"

# This script is intended to be use after CREST calculation.
# The idea is to prune crest_conformers by using invariant RMSD with a 
# defined threshold https://github.com/RMeli/spyrmsd
# After pruning, calculation will be launch using parameters define 
# previously.

# Get cluster specifications
clust = rc.Cluster(fastcar_dir)
# Get keywords for opt calculations
kwds = rc.Keywords(fastcar_dir)
# Get information from parameters.tmp file
paramfile = rp.Parameters('../parameters.tmp')
paramfile.getrmsd()
paramfile.getopt()
paramfile.getwoc()

# Check if CREST terminated normally
crestgood,noreftopo = fm.crest_normterm('crestAnalysis.txt',paramfile.ref_freq)

# Write informations in the summary file
now = datetime.datetime.now()
date = now.strftime("%d/%m/%y %H:%M")
with open(f'../{sumfilename}','a') as sumfile:
    sumfile.write(f'Calling engine.py {date}\n')

if crestgood:
    with open(f'../{sumfilename}','a') as sumfile:
        sumfile.write('CREST terminated normally\n')
    if paramfile.loopn:
        confo_fn = f'crest_conformers_{paramfile.loopn}.xyz'
    else:
        confo_fn = 'crest_conformers.xyz'
    if paramfile.ref_freq:
        shutil.copy('crest_ensemble.xyz',f'../TS-CREST/{confo_fn}') 
        os.chdir('../TS-CREST')
    else:
        shutil.copy('crest_conformers.xyz',f'../Geo-CREST/{confo_fn}') 
        os.chdir('../Geo-CREST')

    # Dict to store geometries of conformers which will be reopt
    crest_uconfos = dict()
    rmsd_thresh = float(paramfile.RMSD_crest)
    reopts = list()
    # Structures from this loop
    comp_mols = io.loadallmols(confo_fn)
    # These ones won't get stripped of their hydrogen atoms
    geo_mols = io.loadallmols(confo_fn)
    nstrucs = len(comp_mols)
    with open(f'../{sumfilename}','a') as sumfile:
        sumfile.write(f'{nstrucs} conformers found by CREST\n')
    # If it is already in a loop it will compare results of new CREST with
    # the results of the previous CREST and won't reoptimize structure close
    # from the ones that have already been optimized
    if paramfile.loopn:
        # Structures from previous loops
        ref_mols = io.loadallmols(uconfilename)
        # RMSD comparison of loop structures with previous structures
        sim = fm.geocomp(comp_mols,ref_mols,rmsd_thresh)
        # Keep only new structures
        ncomp_mols = [cm for c,cm in enumerate(comp_mols) if c not in sim]
        ngeo_mols = [cm for c,cm in enumerate(geo_mols) if c not in sim]
        # Keep in mind the numbering
        sn2 = [c for c,cm in enumerate(comp_mols) if c not in sim]
    else:
        ncomp_mols = comp_mols
        ngeo_mols = geo_mols
        sn2 = range(0,len(comp_mols))
    if len(ncomp_mols) > 0:
        # RMSD comparison between structures from this loop
        sim2 = fm.geocomp_onefile(ncomp_mols,rmsd_thresh)
    else:
        sim2 = []
    # For unique structures, launch optimisation
    for c,cm in enumerate(ngeo_mols):
        if c in sim2:
            continue
        geo = cm.coordinates
        if paramfile.loopn:
            jobname = f"{paramfile.calc_name}-{sn2[c]+1}_loop{paramfile.loopn}"
        else:
            jobname = f"{paramfile.calc_name}-{sn2[c]+1}"
        inpname = f'{jobname}.inp'
        if paramfile.ref_freq:
            if paramfile.guidopt and paramfile.ref_freq:
                calctype = 'opt-ts-guided'
                kw = (kwds.keywords[paramfile.software][calctype],calctype)
                intcoo = list()
                if paramfile.bonds:
                    intcoo.append(paramfile.bonds)
                if paramfile.angles:
                    intcoo.append(paramfile.angles)
                if paramfile.dihedrals:
                    intcoo.append(paramfile.dihedrals)
                fm.writeinp(inpname,kw,geo,cm.atomicnums,paramfile,
                            intcoo=intcoo)
            else:
                calctype = 'opt-ts'
                kw = (kwds.keywords[paramfile.software][calctype],calctype)
                fm.writeinp(inpname,kw,geo,cm.atomicnums,paramfile)
        else:
            calctype = 'opt'
            kw = (kwds.keywords[paramfile.software][calctype],calctype)
            fm.writeinp(inpname,kw,geo,cm.atomicnums,paramfile)
        crest_uconfos[inpname] = (cm.atomicnums,geo)
        if paramfile.loopn:
            reopts.append(f'crestconf{sn2[c]+1}_loop{paramfile.loopn}')
        else:
            reopts.append(f'crestconf{sn2[c]+1}')
        fm.writecalcsubscript(paramfile,jobname,clust)
        if not paramfile.woc:
            os.system(f"{clust.subco} {jobname}.sub")
    if paramfile.loopn:
        pruname = f'pruningCREST_{paramfile.loopn}.log'
    else:
        pruname = f'pruningCREST.log'

    # Write detailled informations about results of pruning in pruning file and 
    # some informations in summary file
    with open(pruname, 'w') as prunf, open(f'../{sumfilename}','a') as sumfile:
        sumfile.write(f'RMSD threshold: {rmsd_thresh}\n')
        sumfile.write(f'{len(reopts)} new conformers\n')
        emp = ''
        heading = 'Results of the pruning' 
        prunf.write(f'{emp:-^80}\n')
        prunf.write(f'{heading:-^80}\n')
        prunf.write(f'{emp:-^80}\n\n')
        prunf.write(f'RMSD threshold: {rmsd_thresh}\n')
        if paramfile.loopn:
            prunf.write(f'{nstrucs} conformers found, {len(reopts)} '
                        'reoptimisations submitted\n\n')
            heading = 'Structures reoptimised' 
        else:
            prunf.write(f'{nstrucs} conformers found, {len(reopts)} '
                        'optimisations submitted\n\n')
            heading = 'Structures optimised' 
        prunf.write(f'\n{heading:-^80}\n')
        for ro in reopts:
            prunf.write(f'{ro}\n')
        if paramfile.loopn:
            if len(sim) > 0:
                sumfile.write(f'{len(sim)} identical to a structure from a '
                               'previous loop\n')
                heading = 'Structure identical to one in previous step' 
                prunf.write(f'\n{heading:-^80}\n')
                with open(uconfilename, 'r') as uconfile:
                    ufclines = uconfile.readlines()
                natoms = int(ufclines[0])
                for s in sim:
                    prevn = sim[s][0]
                    prev_name = ufclines[(natoms+2)*prevn+1].strip('\n')\
                                                            .split('.')[0]
                    prunf.write(f'{s+1}: {prev_name} {sim[s][1]}\n')
        if len(sim2) > 0:
            sumfile.write(f'{len(sim2)} identical to another structure from '
                            'the same loop\n')
            heading = 'Structure identical to one from the same loop' 
            prunf.write(f'\n{heading:-^80}\n')
            for sn in sim2:
               prunf.write(f'{sn2[sn]+1}: {sn2[sim2[sn][0]]+1} {sim2[sn][1]}\n')

    if paramfile.loopn:
        if len(reopts) == 0:
            # For loop calculations, if no new structures, stop the loops and
            # go to abs.py
            with open(f'../{sumfilename}','a') as sumfile:
                sumfile.write('No new structure found by CREST in check '
                              'calculation\n')
                sumfile.write(f'{emp:-^80}\n{emp:-^80}\n\n')
            os.chdir('../Final/')
            scriptname = fm.writesubscript('abs',paramfile,clust)
            os.system(f"{clust.subco} {scriptname}")   
            sys.exit(0)
    # Write opt parameters in summary file
    with open(f'../{sumfilename}', 'a') as sumfile:
        emp = ''
        sumfile.write(f'{emp:-^80}\n{emp:-^80}\n\n')
        heading = 'Optimisations'
        sumfile.write(f'{emp:-^80}\n')
        sumfile.write(f'{heading:-^80}\n\n')
        heading = 'Opt parameters'
        sumfile.write(f'{heading:-^80}\n')
        sumfile.write(f'Software: {paramfile.software} \n')
        sumfile.write(f'Method: {paramfile.method}\n')
        sumfile.write(f'Basis set: {paramfile.basis_set}\n')
        if paramfile.dispersion != '':
            sumfile.write(f'Dispersion: {paramfile.dispersion}\n')
        if paramfile.solvent != 'none':
            sumfile.write(f'Solvent: {paramfile.solvent}\n')
        sumfile.write(f'\n')
        sumfile.write(f'{len(reopts)} gaussian optimization submitted\n')
    if paramfile.loopn:
        wm = 'a'
    else:
        wm = 'w'
    # Write or update file with all unique CREST conformers
    fm.writegeolist(uconfilename,crest_uconfos,wm)
    # Call transmission script that will start when all reoptimisations are done
    if paramfile.woc:
        dep = False
    else:
        dep = True
    scriptname = fm.writesubscript('transmission',paramfile,clust,
                                   dependency=dep)
    os.system(f'{clust.subco} {scriptname}')

elif noreftopo:
    # Check if CREST did not proceed due to topology error and relaunch it 
    # using noreftopo
    with open(f'../{sumfilename}','a') as sumfile:
        sumfile.write('Topology error in CREST, resubmitting calculation with'
                      '--noreftopo\n')
    if loop:
        crestjobname = f'crest{paramfile.loopn}_{startoutput.base_name}'
    else:
        crestjobname= f'crest_{startoutput.base_name}'
    with open(crestjobname, 'r') as file:
        filedata = file.read()
    if paramfile.ref_freq == '':
        filedata = filedata.replace('> crestAnalysis.txt', '--noreftopo > '
                                    'crestAnalysis.txt')
    else :
        filedata = filedata.replace('-cinp constraints.inp --subrmsd >'
                                    ' crestAnalysis.txt', 
                                    '--noreftopo -cinp constraints.inp '
                                    '--subrmsd > crestAnalysis.txt')
    with open(crestjobname, 'w') as file:
        file.write(filedata)
    os.system(f"{clust.subco} {crestjobname}")
    with open(f'../{sumfilename}', 'a') as out:
        out.write('CREST restarted with noreftopo\n')

else:
   # If CREST failed, stop 
    with open(f'../{sumfilename}','a') as sumfile:
        emp = ''
        sumfile.write('CREST failed\n')
        sumfile.write(f'{emp:-^80}\n{emp:-^80}\n\n')


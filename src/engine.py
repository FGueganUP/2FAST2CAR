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
# threshold set to 0.5. https://github.com/RMeli/spyrmsd
# After pruning, calculation will be launch using DFT parameters define 
# previously.

# Get cluster specifications
clust = rc.Cluster(fastcar_dir)
# Get keywords
kwds = rc.Keywords(fastcar_dir)
# Get information from parameters_mod.txt file
paramfile = rp.Parameters('../parameters.tmp')

# Check if CREST terminated normally
crestgood,noreftopo = fm.crest_normterm('CrestAnalysis.txt',paramfile.ref_freq)

now = datetime.datetime.now()
with open(f'../{sumfilename}','a') as sumfile:
    sumfile.write(f'Calling engine.py at {now.time()}\n')

if crestgood:
    with open(f'../{sumfilename}','a') as sumfile:
        sumfile.write('CREST terminated normally\n')
    if paramfile.check_st:
        confo_fn = f'crest_conformers_{paramfile.check_st}.xyz'
    else:
        confo_fn = 'crest_conformers.xyz'
    if paramfile.ref_freq:
        #Create directory for TS calculations
        if not paramfile.woc:
            if os.path.isdir('../TS-CREST'):
                if not paramfile.check_st:
                    shutil.rmtree('../TS-CREST')
                    os.mkdir('../TS-CREST')
            else:
                os.mkdir('../TS-CREST')
        shutil.copy('crest_ensemble.xyz',f'../TS-CREST/{confo_fn}') 
        os.chdir('../TS-CREST')
    else:
         #Create directory for geometry calculations  
        if not paramfile.woc:
            if os.path.isdir('../Geo-CREST'):
                if not paramfile.check_st:
                    shutil.rmtree('../Geo-CREST')
                    os.mkdir('../Geo-CREST')
            else:
                os.mkdir('../Geo-CREST')
        shutil.copy('crest_conformers.xyz',f'../Geo-CREST/{confo_fn}') 
        os.chdir('../Geo-CREST')
    # Start pruning
    #if paramfile.check_st:
    #    splitname = f'xyzfile{paramfile.check_st}num'
    #else:
    #    splitname = 'xyzfilenum'
    #fm.xyzspliter(confo_fn,splitname)

    #dirs = os.listdir()
    #xyzlist = [ele for ele in dirs if ".xyz" in ele and "num" in ele]
    #ordxyzlist = sorted(xyzlist)
    #nstrucs = len(xyzlist)
    

    #n = 1
    #listn = ordxyzlist[:]

    # Dict to store geometries of conformers which will be reopt
    crest_uconfos = dict()
    # If it is a check calculation it will compare results of new CREST with
    # the results of the previous CREST and won't reoptimize structure close
    # from the ones that have already been optimized
    rmsd_thresh = float(paramfile.RMSD_crest)
    reopts = list()
    comp_mols = io.loadallmols(confo_fn)
    geo_mols = io.loadallmols(confo_fn)
    nstrucs = len(comp_mols)
    with open(f'../{sumfilename}','a') as sumfile:
        sumfile.write(f'{nstrucs} conformers found by CREST\n')
    if paramfile.check_st:
        ref_mols = io.loadallmols(uconfilename)
        sim = fm.geocomp(comp_mols,ref_mols,rmsd_thresh)
        ncomp_mols = [cm for c,cm in enumerate(comp_mols) if c not in sim]
        ngeo_mols = [cm for c,cm in enumerate(geo_mols) if c not in sim]
        sn2 = [c for c,cm in enumerate(comp_mols) if c not in sim]
    else:
        ncomp_mols = comp_mols
        ngeo_mols = geo_mols
        sn2 = range(0,len(comp_mols))
    if len(ncomp_mols) > 0:
        sim2 = fm.geocomp_onefile(ncomp_mols,rmsd_thresh)
    else:
        sim2 = []
    for c,cm in enumerate(ngeo_mols):
        if c in sim2:
            continue
        geo = cm.coordinates
        if paramfile.check_st:
            jobname = f"{paramfile.base_name_1}_{paramfile.experience_number}"\
                      f"-{sn2[c]+1}_loop{paramfile.check_st}"
        else:
            jobname = f"{paramfile.base_name_1}_{paramfile.experience_number}"\
                      f"-{sn2[c]+1}"
        inpname = f'{jobname}.inp'
        soft = 'gaussian'
        if paramfile.ref_freq:
            calctype = 'opt-ts'
        else:
            calctype = 'opt'
        kw = kwds.keywords[soft][calctype]
        fm.writeg16inp(inpname,paramfile,12,geo,cm.atomicnums,kw)
        crest_uconfos[inpname] = (cm.atomicnums,geo)
        if paramfile.check_st:
            reopts.append(f'crestconf{sn2[c]+1}_loop{paramfile.check_st}')
        else:
            reopts.append(f'crestconf{sn2[c]+1}')
        # Create sub file for the unique conformers previously found
        fm.writeg16subscript(paramfile,jobname,clust)
        if not paramfile.woc:
            os.system(f"sbatch {jobname}g16.sub")
    if paramfile.check_st:
        pruname = f'pruningCREST_{paramfile.check_st}.log'
    else:
        pruname = f'pruningCREST.log'
    with open(pruname, 'w') as prunf, open(f'../{sumfilename}','a') as sumfile:
        sumfile.write(f'RMSD threshold: {rmsd_thresh}\n')
        sumfile.write(f'{len(reopts)} new conformers\n')
        emp = ''
        heading = 'Results of the pruning' 
        prunf.write(f'{emp:-^80}\n')
        prunf.write(f'{heading:-^80}\n')
        prunf.write(f'{emp:-^80}\n\n')
        prunf.write(f'RMSD threshold: {rmsd_thresh}\n')
        if paramfile.check_st:
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
        if paramfile.check_st:
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
    #while len(listn) > 0:
    #    if paramfile.check_st:
    #        prevs = (sorted(glob.glob(f'xyzprev*')))
    #        result = subprocess.check_output(f'python -m spyrmsd -m {listn[0]}'
    #                                         ' xyzprev*', shell=True)
    #        rmsd_values = result.split()
    #        rmsd_values = [item.decode() for item in rmsd_values]
    #        recalc = True
    #        for r,rmsd in enumerate(rmsd_values):
    #            if float(rmsd) < float(rmsd_thresh):
    #                recalc = False
    #                already_in_prev[listn[0]] = (prevs[r],rmsd)
    #                break
    #        if recalc:
    #            recalcs += 1
    #        else:
    #            os.remove(listn[0])
    #            listn = listn[1:]
    #            continue
    #    match = re.findall(r'\d+', listn[0])[-1]
    #    # Create gaussian input for unique conformers found
    #    if paramfile.check_st:
    #        jobname = f"{paramfile.base_name_1}_{paramfile.experience_number}"\
    #                  f"-{match}_loop{paramfile.check_st}"
    #    else:
    #        jobname = f"{paramfile.base_name_1}_{paramfile.experience_number}"\
    #                  f"-{match}"
    #    inpname = f'{jobname}.inp'
    #    with open(listn[0], 'r') as fp:
    #        # Extract coordinates from xyz file
    #        geo = fp.read().splitlines(True)[2:]
    #    # Create input for Gaussian calc
    #    soft = 'gaussian'
    #    if paramfile.ref_freq:
    #        calctype = 'opt-ts'
    #    else:
    #        calctype = 'opt'
    #    kw = kwds.keywords[soft][calctype]
    #    fm.writeg16inp(inpname,paramfile,12,geo,kw)
    #    crest_uconfos[inpname] = geo
    #    # Remove unique geom and duplicate and continue the pruning
    #    removals,remov_rmsds=fm.xyzlistcleaner(listn,rmsd_thresh)
    #    #print(f"{removals} removed")
    #    for f,file_to_remove in enumerate(removals):
    #        os.remove(file_to_remove)
    #        id_to_another[file_to_remove] = (listn[0],remov_rmsds[f])
    #    listn = [item for item in listn if item not in removals]
    #    reopts.append(listn[0])
    #    os.remove(listn[0])
    #    listn = listn[1:]

    #    # Create sub file for the unique conformers previously found
    #    fm.writeg16subscript(paramfile,jobname,clust)
    #    if not paramfile.woc:
    #        os.system(f"sbatch {jobname}g16.sub")
    #    n += 1

    #if paramfile.check_st:
    #    pruname = f'pruningCREST_{paramfile.check_st}.log'
    #else:
    #    pruname = f'pruningCREST.log'
    #with open(pruname, 'w') as prunf, open(f'../{sumfilename}','a') as sumfile:
    #    sumfile.write(f'RMSD threshold: {rmsd_thresh}\n')
    #    sumfile.write(f'{len(reopts)} new conformers\n')
    #    emp = ''
    #    heading = 'Results of the pruning' 
    #    prunf.write(f'{emp:-^80}\n')
    #    prunf.write(f'{heading:-^80}\n')
    #    prunf.write(f'{emp:-^80}\n\n')
    #    prunf.write(f'RMSD threshold: {rmsd_thresh}\n')
    #    if paramfile.check_st:
    #        prunf.write(f'{nstrucs} conformers found, {len(reopts)} '
    #                    'reoptimisations submitted\n\n')
    #        heading = 'Structures reoptimised' 
    #    else:
    #        prunf.write(f'{nstrucs} conformers found, {len(reopts)} '
    #                    'optimisations submitted\n\n')
    #        heading = 'Structures optimised' 
    #    prunf.write(f'\n{heading:-^80}\n')
    #    for ro in reopts:
    #        prunf.write(f'{ro}\n')
    #    if len(already_in_prev) > 0:
    #        sumfile.write(f'{len(already_in_prev)} identical to a structure '
    #                      'from a previous loop\n')
    #        heading = 'Structure identical to one in previous step' 
    #        prunf.write(f'\n{heading:-^80}\n')
    #        for aip in already_in_prev:
    #            prunf.write(f'{aip}: {already_in_prev[aip]}\n')
    #    if len(id_to_another) > 0:
    #        sumfile.write(f'{len(id_to_another)} identical to another '
    #                      'structure from the same loop\n')
    #        heading = 'Structure identical to one already reoptimised' 
    #        prunf.write(f'\n{heading:-^80}\n')
    #        for ita in id_to_another:
    #            prunf.write(f'{ita}: {id_to_another[ita]}\n')

    if paramfile.check_st:
        if len(reopts) == 0:
            # For check calculations, if no new structures, stop the process
            with open(f'../{sumfilename}','a') as sumfile:
                sumfile.write('No new structure found by CREST in check '
                              'calculation\n')
            os.chdir('../Final/')
            fm.writeabsscript(paramfile)
            os.system(f"sbatch script_abs.sub")   
            sys.exit(0)
    with open(f'../{sumfilename}', 'a') as sumfile:
        emp = ''
        heading = 'DFT parameters'
        sumfile.write(f'\n{heading:-^80}\n')
        sumfile.write(f'Functional: {paramfile.functional}\n')
        sumfile.write(f'Basis set: {paramfile.base_name_1}\n')
        sumfile.write(f'Dispersion: {paramfile.dispersion}\n')
        sumfile.write(f'Solvent: {paramfile.solvent}\n')
        sumfile.write(f'{emp:-^80}\n\n')
        sumfile.write(f'{len(reopts)} gaussian optimization submitted\n')
    if paramfile.check_st:
        wm = 'a'
    else:
        wm = 'w'
    fm.writegeolist(uconfilename,crest_uconfos,wm)
    # Submit the following script that will analyse previously DFT calculations
    fm.writecheckscript(paramfile)
    os.system(f'sbatch script_check.sub')

elif noreftopo:
    # Check if CREST did not proceed due to topology error and relaunch it 
    # using noreftopo
    with open(f'../{sumfilename}','a') as sumfile:
        sumfile.write('Topology error in CREST, resubmitting calculation with'
                      '--noreftopo\n')
    with open('script_CREST.sub', 'r') as file:
        filedata = file.read()
    if paramfile.ref_freq == '':
        filedata = filedata.replace('> CrestAnalysis.txt', '--noreftopo > '
                                    'CrestAnalysis.txt')
    else :
        filedata = filedata.replace('-cinp constraints.inp --subrmsd >'
                                    ' CrestAnalysis.txt', 
                                    '--noreftopo -cinp constraints.inp '
                                    '--subrmsd > CrestAnalysis.txt')
    with open('script_CREST.sub', 'w') as file:
        file.write(filedata)
    os.system(f"sbatch script_CREST.sub")
    with open(f'../{sumfilename}', 'a') as out:
        out.write('CREST restarted with noreftopo\n')

else:
    with open(f'../{sumfilename}','a') as sumfile:
        sumfile.write('CREST failed\n')


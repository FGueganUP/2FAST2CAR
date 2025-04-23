#!/usr/bin/env python3
import sys
import re
import glob
import os
import datetime
import shutil
import subprocess
import math
import numpy as np
import scipy as sp
import spyrmsd
from spyrmsd import io, rmsd
fastcar_dir = "/home/hgeindre/ProjectFASTCAR/FASTCAR2"
sys.path.append(fastcar_dir)
import readparam as rp
import readoutput as ro
import filemanager as fm
import readconfig as rc
import misc

# This script is used after optimisations, it will check for failed 
# calculations, duplicates, wrong TS...

prevfilename = 'prevOPT.xyz'
loopfilename = 'loopOPT.xyz'
sumfilename = 'summary.log'

# Get cluster specifications
clust = rc.Cluster(fastcar_dir)
# Get information from parameters.tmp file
paramfile = rp.Parameters('../parameters.tmp')
paramfile.getopt()
paramfile.getwoc()
paramfile.getloop()
paramfile.getrmsd()

now = datetime.datetime.now()
date = now.strftime("%d/%m/%y %H:%M")
with open(f'../{sumfilename}','a') as sumfile:
    sumfile.write(f'Calling transmission {date}\n')

# Case of TS
if paramfile.ref_freq:
    uniqueTS = []
    duplicateTS = []
    garbageTS = []
    noTS = []
    failedTS = []
    
    if paramfile.loopn:
        loopfiles =  glob.glob(f'{paramfile.calc_name}*_'
                               f'loop{paramfile.loopn}.out')
        compTSname = f'compTS_{paramfile.loopn}.log'
        file_name = glob.glob('../Final/*_lowest.out')[0]
        startoutput = ro.Output(f'../Final/{file_name}')
    else:
        loopfiles = glob.glob(f'{paramfile.calc_name}*.out')
        compTSname = f'compTS.log'
        startoutput = ro.Output(f'../{paramfile.calc_name}.out')
    goodfiles = list()
    #Create file to write information on TS screening
    if paramfile.tsc == 'scalprod':
        compTS =  open(compTSname,'w')
        emp = ''
        heading = 'Scalar product of reaction normal modes' 
        compTS.write(f'{emp:-^80}\n')
        compTS.write(f'{heading:-^80}\n')
        compTS.write(f'{emp:-^80}\n\n')

    #Check result of reoptimisations
    for filename in loopfiles:
        outputfile = ro.Output(filename)
        # Check for failed calculations
        if not outputfile.normterm:
            failedTS.append((filename, 'Calculation failed'))
            # Change the name of the failed TS
            if not 'failed' in filename:
                new_name = filename.split('.')[0] + '_failed.' + \
                           filename.split('.')[1]
                shutil.move(filename,new_name)
            continue
        x = outputfile.freqs[0]
        # Check for failed freq calculations
        if math.isnan(x):
            # Change the name of the failed TS
            if not 'failed' in filename:
                new_name = filename.split('.')[0] + '_failed.' + \
                           filename.split('.')[1]
                shutil.move(filename,new_name)
            else:
                new_name = filename
            failedTS.append((new_name, 'Freq calculation failed'))
            continue
        # Check if TS
        if outputfile.ts:
            # Check if TS describes the same reaction as the start file
            if paramfile.tsc == 'activats':
                isid = misc.checknv(outputfile.normvec,paramfile.activ_ats)
            elif paramfile.tsc == 'scalprod':
                rot,rssd=sp.spatial.transform.Rotation.align_vectors\
                         (startoutput.coordinates,outputfile.coordinates)
                rot_mat = rot.as_matrix()
                anormvec = list()
                for v in outputfile.normvec[0]:
                    new_vec = rot_mat.dot(v)
                    anormvec.append(new_vec)
                isid1,scalp1 = misc.compnv(startoutput.normvec[0],anormvec)
                isid2,scalp2 = misc.compnv(startoutput.normvec[0],
                                           outputfile.normvec[0])
                if scalp1 > scalp2:
                    scalp = scalp1
                    isid = isid1
                else:
                    scalp = scalp2
                    isid = isid2
                compTS.write(f'{filename}: {scalp1} {scalp2}\n')
            if isid:
                goodfiles.append(filename)
            else:
                # Change the name of the garbage outputfiles
                if not 'garbage' in filename:
                    new_name = filename.split('.')[0] + '_garbage.' + \
                               filename.split('.')[1]
                    shutil.move(filename,new_name)
                else:
                    new_name = filename
                garbageTS.append((new_name, float(outputfile.last_SCF), x))
        else:
            # Change the name of the outputfiles that are not TS
            if not 'noTS' in filename:
                new_name = filename.split('.')[0] + '_noTS.' + \
                           filename.split('.')[1]
            else:
                new_name = filename
            shutil.move(filename,new_name)
            noTS.append((new_name, float(outputfile.last_SCF), x))

    rmsd_thresh = paramfile.RMSD_opt
    id_to_another = dict()
    newstrucs = list()
    duplicates = list()
    # Check if opt structures are identical to ones from previous loops
    if paramfile.loopn:
        prevfiles = glob.glob(f'{paramfile.calc_name}*.out')
        prevfiles = [fn for fn in prevfiles if not f'loop{paramfile.loopn}'\
                     in fn and not 'duplicate' in fn and not 'failed' in fn \
                     and not 'garbage' in fn and not 'noTS' in fn]
        fm.outstoxyzfile(prevfiles,prevfilename)
        # Opt structures from previous loops
        prev_mols = io.loadallmols(prevfilename)
        fm.outstoxyzfile(goodfiles,loopfilename)
        # Opt structures from this loop
        loop_mols = io.loadallmols(loopfilename)
        # These ones won't have their hydrogen stripped
        geo_mols = io.loadallmols(loopfilename)
        sim = fm.geocomp(loop_mols,prev_mols,rmsd_thresh)
        nloop_mols = [cm for c,cm in enumerate(loop_mols) if c not in sim]
        ngeo_mols = [cm for c,cm in enumerate(geo_mols) if c not in sim]
        # Keep numbering in mind
        sn2 = [c for c,cm in enumerate(loop_mols) if c not in sim]
    else:
        fm.outstoxyzfile(goodfiles,loopfilename)
        nloop_mols = io.loadallmols(loopfilename)
        sn2 = range(0,len(nloop_mols))
    if len(nloop_mols) > 0:
        # Check for identical structures in the ones from this loop
        sim2 = fm.geocomp_onefile(nloop_mols,rmsd_thresh)
    else:
        sim2 = []
    if paramfile.loopn:
        # Structure identical to one from a previous loop
        for s in sim:
            filename = goodfiles[s]
            duplicates.append(filename)
            # Change the name of the duplicates
            if not 'duplicate' in filename:
                outname = filename.split('.')[0] + '.out'
                new_name = outname.split('.')[0] + '_duplicate.out'
                shutil.move(outname,new_name)
            else:
                new_name = filename
            id_to_another[new_name] = (prevfiles[sim[s][0]], sim[s][1])
    # Structure identical to one from the same loop
    for s in sim2:
        filename = goodfiles[sn2[s]]
        duplicates.append(filename)
        # Change the name of the duplicates
        if not 'duplicate' in filename:
            outname = filename.split('.')[0] + '.out'
            new_name = outname.split('.')[0] + '_duplicate.out'
            shutil.move(outname,new_name)
        else:
            new_name = filename
        id_to_another[new_name] = (goodfiles[sn2[sim2[s][0]]], sim2[s][1])

    for fn in goodfiles:
        if fn not in duplicates:
            newstrucs.append(fn)
        
    # Good TSs
    for filename in newstrucs:
        outputfile = ro.Output(filename)
        uniqueTS.append((filename, float(outputfile.last_SCF),
                         outputfile.freqs[0]))
            
    compname = f'{paramfile.calc_name}*.out'

    # Write results of opt pruning in pruning file
    if paramfile.loopn:
        pruname = f'pruningOPT_{paramfile.loopn}.log'
    else:
        pruname = f'pruningOPT.log'
    with open(pruname,'w') as prunf:
        emp = ''
        heading = 'Results of the pruning' 
        prunf.write(f'{emp:-^80}\n')
        prunf.write(f'{heading:-^80}\n')
        prunf.write(f'{emp:-^80}\n\n')
        prunf.write(f'RMSD threshold: {rmsd_thresh}\n')
        prunf.write(f'{len(loopfiles)} optimisation done, {len(newstrucs)} new'
                    ' structures\n\n')
        heading = 'New structures' 
        for ns in newstrucs:
            prunf.write(f'{ns}\n')
        if len(id_to_another) > 0:
            heading = 'Structure identical to another one' 
            prunf.write(f'\n{heading:-^80}\n')
            for ita in id_to_another:
                prunf.write(f'{ita}: {id_to_another[ita]}\n')
    for filename in id_to_another:
        outputfile = ro.Output(filename)
        duplicateTS.append((filename, outputfile.last_SCF, outputfile.freqs[0],
                         id_to_another[filename][0],id_to_another[filename][1]))

    #Sort all dict containing the different results of optimisation by energy
    uniqueTS.sort(key=lambda item: item[1])
    duplicateTS.sort(key=lambda item: item[1])
    garbageTS.sort(key=lambda item: item[1])
    noTS.sort(key=lambda item: item[1])
    failedTS.sort(key=lambda item: item[0])
    if paramfile.loopn:
        wm = 'a'
    else:
        wm = 'w'
    # Update or create structures.log file with informations about reoptimised 
    # structures
    fm.writelogTS(uniqueTS,duplicateTS,garbageTS,noTS,failedTS,wm)

    # Write informations about results of opt calculation in summary file
    with open(f'../{sumfilename}', 'a') as sumfile:
        sumfile.write('Informations about optimised TS written in '
                      'structures.log\n\n')
        sumfile.write(f'RMSD threshold: {rmsd_thresh}\n')
        sumfile.write(f'{len(uniqueTS)} new TS(s)\n') 
        sumfile.write(f'{len(duplicateTS)} duplicate(s)\n') 
        sumfile.write(f'{len(garbageTS)} TS(s) describing another reaction\n') 
        sumfile.write(f'{len(noTS)} that are not TS\n') 
        sumfile.write(f'{len(failedTS)} calculation(s) failed\n') 

    # Decide what to do next, continue loops, go to abs.py or stop calculation
    # If no loops, go to abs.py to finalize the process
    if paramfile.loop == 'none':
        if not paramfile.woc:
            if os.path.isdir('../Final/'):
                shutil.rmtree('../Final/')
            os.mkdir('../Final/')
        if len(uniqueTS) > 0:
            lowestname = uniqueTS[0][0].split('.')[0] + '_lowest.' +\
                         uniqueTS[0][0].split('.')[1]
            shutil.copy(uniqueTS[0][0], f'../Final/{lowestname}')
        else:
            with open(f'../{sumfilename}','a') as sumfile:
                sumfile.write(f"No conformer found, ending calculation\n")
                sumfile.write(f'{emp:-^80}\n{emp:-^80}\n\n')
                sys.exit(0)
        os.chdir('../Final/')
        scriptname = fm.writesubscript('abs',paramfile,clust)
        os.system(f"{clust.subco} {scriptname}")   
    # If in first loop, start a new one
    elif paramfile.loopn == 0:
        if not paramfile.woc:
            if os.path.isdir('../Final/'):
                shutil.rmtree('../Final/')
            os.mkdir('../Final/')
        if len(uniqueTS) > 0:
            lowestname = uniqueTS[0][0].split('.')[0] + '_lowest.' +\
                         uniqueTS[0][0].split('.')[1]
            shutil.copy(uniqueTS[0][0], f'../Final/{lowestname}')
        else:
            with open(f'../{sumfilename}','a') as sumfile:
                sumfile.write(f"No conformer found, ending calculation\n")
                sumfile.write(f'{emp:-^80}\n{emp:-^80}\n\n')
                sys.exit(0)
        os.chdir('../')
        scriptname = fm.writesubscript('sparkplug',paramfile,clust)
        os.system(f"{clust.subco} {scriptname}")   
    # If already looping, check if new opt structures were found and if number 
    # of max loops has not been exceeded, continue looping if so. If not, call
    # abs.py
    else:
        if len(uniqueTS) > 0:
            plname = glob.glob('../Final/*_lowest.out')[0]
            prevlowest = ro.Output(plname)
            if uniqueTS[0][2] < prevlowest.last_SCF:
                lowestname = uniqueTS[0][0].split('.')[0] + '_lowest.' +\
                             uniqueTS[0][0].split('.')[1]
                shutil.copy(f'{uniqueTS[0][0]}', f'../Final/{lowestname}')
                os.remove(plname)
                with open(f'../{sumfilename}','a') as sumfile:
                    sumfile.write(f"New lowest energy structure: {lowestname}\n")
                    sumfile.write(f'{emp:-^80}\n{emp:-^80}\n\n')
            else:
                with open(f'../{sumfilename}','a') as sumfile:
                    sumfile.write(f"No new lowest energy structure\n")
                    sumfile.write(f'{emp:-^80}\n{emp:-^80}\n\n')
            if paramfile.loopn == paramfile.loop:
                os.chdir('../Final/')
                scriptname = fm.writesubscript('abs',paramfile,clust)
                os.system(f"{clust.subco} {scriptname}")   
            else:
                os.chdir('../')
                scriptname = fm.writesubscript('sparkplug',paramfile,clust)
                os.system(f"{clust.subco} {scriptname}")   
        else:
            with open(f'../{sumfilename}','a') as sumfile:
                sumfile.write("No new structure\n")
                sumfile.write(f'{emp:-^80}\n{emp:-^80}\n\n')
            os.chdir('../Final/')
            scriptname = fm.writesubscript('abs',paramfile,clust)
            os.system(f"{clust.subco} {scriptname}")   

# Case of minimum
else:
    uniqueMin = []
    duplicateMin = []
    noMin = []
    failedMin = []
    if paramfile.loopn:
        loopfiles =  glob.glob(f'{paramfile.calc_name}*_'
                               f'loop{paramfile.loopn}.out')
    else:
        loopfiles = glob.glob(f'{paramfile.calc_name}*.out')

    if paramfile.loopn:
        file_name = glob.glob('../Final/*_lowest.out')[0]
        startoutput = ro.Output(f'../Final/{file_name}')
    else:
        startoutput = ro.Output(f'../{paramfile.calc_name}.out')
    goodfiles = list()
    for filename in loopfiles:
        outputfile = ro.Output(filename)
        if not outputfile.normterm:
            failedMin.append((filename, 'Calculation failed'))
            # Change the name of the failed calculation
            if not 'failed' in filename:
                new_name = filename.split('.')[0] + '_failed.' + \
                           filename.split('.')[1]
                shutil.move(filename,new_name)
            continue
        x = outputfile.freqs[0]
        if math.isnan(x):
            # Change the name of the failed calculation
            if not 'failed' in filename:
                new_name = filename.split('.')[0] + '_failed.' + \
                           filename.split('.')[1]
                shutil.move(filename,new_name)
            else:
                new_name = filename
            failedMin.append((new_name, 'Freq calculation failed'))
            continue
        if outputfile.ts:
            if not 'garbage' in filename:
                new_name = filename.split('.')[0] + '_garbage.' + \
                           filename.split('.')[1]
                shutil.move(filename,new_name)
            else:
                new_name = filename
            noMin.append((new_name, outputfile.last_SCF, x))
            continue
        goodfiles.append(filename)
    rmsd_thresh = 0.1
    id_to_another = dict()
    newstrucs = list()
    duplicates = list()
    # Check if opt structures are identical to ones from previous loops
    if paramfile.loopn:
        prevfiles = glob.glob(f'{paramfile.calc_name}*.out')
        prevfiles = [fn for fn in prevfiles if not f'loop{paramfile.loopn}'\
                     in fn and not 'duplicate' in fn and not 'failed' in fn \
                     and not 'garbage' in fn]
        fm.outstoxyzfile(prevfiles,prevfilename)
        # Opt structures from previous loops
        prev_mols = io.loadallmols(prevfilename)
        fm.outstoxyzfile(goodfiles,loopfilename)
        # Opt structures from this loop
        loop_mols = io.loadallmols(loopfilename)
        # These ones won't have their hydrogen stripped
        geo_mols = io.loadallmols(loopfilename)
        sim = fm.geocomp(loop_mols,prev_mols,rmsd_thresh)
        nloop_mols = [cm for c,cm in enumerate(loop_mols) if c not in sim]
        ngeo_mols = [cm for c,cm in enumerate(geo_mols) if c not in sim]
        # Keep numbering in mind
        sn2 = [c for c,cm in enumerate(loop_mols) if c not in sim]
    else:
        fm.outstoxyzfile(goodfiles,loopfilename)
        nloop_mols = io.loadallmols(loopfilename)
        sn2 = range(0,len(nloop_mols))
    if len(nloop_mols) > 0:
        # Check for identical structures in the ones from this loop
        sim2 = fm.geocomp_onefile(nloop_mols,rmsd_thresh)
    else:
        sim2 = []
    if paramfile.loopn:
        # Structure identical to one from a previous loop
        for s in sim:
            filename = goodfiles[s]
            duplicates.append(filename)
            # Change the name of the duplicates
            if not 'duplicate' in filename:
                outname = filename.split('.')[0] + '.out'
                new_name = outname.split('.')[0] + '_duplicate.out'
                shutil.move(outname,new_name)
            else:
                new_name = filename
            id_to_another[new_name] = (prevfiles[sim[s][0]], sim[s][1])
    # Structure identical to one from the same loop
    for s in sim2:
        filename = goodfiles[sn2[s]]
        duplicates.append(filename)
        # Change the name of the duplicates
        if not 'duplicate' in filename:
            outname = filename.split('.')[0] + '.out'
            new_name = outname.split('.')[0] + '_duplicate.out'
            shutil.move(outname,new_name)
        else:
            new_name = filename
        id_to_another[new_name] = (goodfiles[sn2[sim2[s][0]]], sim2[s][1])

    for fn in goodfiles:
        if fn not in duplicates:
            newstrucs.append(fn)
    # Good mins
    for filename in newstrucs:
        outputfile = ro.Output(filename)
        uniqueMin.append((filename, outputfile.last_SCF))
            
    compname = f'{paramfile.calc_name}*.out'

    # Write results of opt pruning in pruning file
    if paramfile.loopn:
        pruname = f'pruningOPT_{paramfile.loopn}.log'
    else:
        pruname = f'pruningOPT.log'
    with open(pruname,'w') as prunf:
        emp = ''
        heading = 'Results of the pruning' 
        prunf.write(f'{emp:-^80}\n')
        prunf.write(f'{heading:-^80}\n')
        prunf.write(f'{emp:-^80}\n\n')
        prunf.write(f'RMSD threshold: {rmsd_thresh}\n')
        prunf.write(f'{len(loopfiles)} optimisation done, {len(newstrucs)} new'
                    ' structures\n\n')
        heading = 'New structures' 
        for ns in newstrucs:
            prunf.write(f'{ns}\n')
        if len(id_to_another) > 0:
            heading = 'Structure identical to another one' 
            prunf.write(f'\n{heading:-^80}\n')
            for ita in id_to_another:
                prunf.write(f'{ita}: {id_to_another[ita]}\n')
    for filename in id_to_another:
        outputfile = ro.Output(filename)
        duplicateMin.append((filename,outputfile.last_SCF,
                           id_to_another[filename][0],
                           id_to_another[filename][1]))

    uniqueMin.sort(key=lambda item: item[1])
    duplicateMin.sort(key=lambda item: item[1])
    failedMin.sort(key=lambda item: item[0])
    noMin.sort(key=lambda item: item[1])
    # Write in log file sorted geom by category
    if paramfile.loopn:
        wm = 'a'
    else:
        wm = 'w'
    fm.writelogmin(uniqueMin,duplicateMin,noMin,failedMin,wm)
    with open(f'../{sumfilename}', 'a') as sumfile:
        sumfile.write('Informations about optimised minima written in '
                      'structures.log\n\n')
        sumfile.write(f'RMSD threshold: {rmsd_thresh}\n')
        sumfile.write(f'{len(uniqueMin)} new minimum(a)\n') 
        sumfile.write(f'{len(duplicateMin)} duplicate(s)\n') 
        sumfile.write(f'{len(noMin)} with a negative frequency\n') 
        sumfile.write(f'{len(failedMin)} optimisation(s) failed\n') 

    # If no loops, go to abs.py to finalize the process
    if paramfile.loop == 'none':
        if not paramfile.woc:
            if os.path.isdir('../Final/'):
                shutil.rmtree('../Final/')
            os.mkdir('../Final/')
        if len(uniqueMin) > 0:
            lowestname = uniqueMin[0][0].split('.')[0] + '_lowest.' +\
                         uniqueMin[0][0].split('.')[1]
            shutil.copy(uniqueMin[0][0], f'../Final/{lowestname}')
        else:
            with open(f'../{sumfilename}','a') as sumfile:
                sumfile.write(f"No conformer found, ending calculation\n")
                sumfile.write(f'{emp:-^80}\n{emp:-^80}\n\n')
                sys.exit(0)
        os.chdir('../Final/')
        scriptname = fm.writesubscript('abs',paramfile,clust)
        os.system(f"{clust.subco} {scriptname}")   
    # If in first loop, start a new one
    elif paramfile.loopn == 0:
        if not paramfile.woc:
            if os.path.isdir('../Final/'):
                shutil.rmtree('../Final/')
            os.mkdir('../Final/')
        if len(uniqueMin) > 0:
            lowestname = uniqueMin[0][0].split('.')[0] + '_lowest.' +\
                         uniqueMin[0][0].split('.')[1]
            shutil.copy(uniqueMin[0][0], f'../Final/{lowestname}')
        else:
            with open(f'../{sumfilename}','a') as sumfile:
                sumfile.write(f"No conformer found, ending calculation\n")
                sys.exit(0)
        os.chdir('../')
        scriptname = fm.writesubscript('sparkplug',paramfile,clust)
        os.system(f"{clust.subco} {scriptname}")   
    # If already looping, check if new opt structures were found and continue
    # looping if so. If not, go to abs.py
    else:
        if len(uniqueMin) > 0:
            plname = glob.glob('../Final/*_lowest.out')[0]
            prevlowest = ro.Output(plname)
            if uniqueMin[0][1] < prevlowest.last_SCF:
                lowestname = uniqueMin[0][0].split('.')[0] + '_lowest.' +\
                             uniqueMin[0][0].split('.')[1]
                shutil.copy(f'{uniqueMin[0][0]}', f'../Final/{lowestname}')
                os.remove(plname)
                with open(f'../{sumfilename}','a') as sumfile:
                    sumfile.write("New lowest energy structure: "
                                 f"{lowestname}\n")
                    sumfile.write(f'{emp:-^80}\n{emp:-^80}\n\n')
            else:
                with open(f'../{sumfilename}','a') as sumfile:
                    sumfile.write(f"No new lowest energy structure\n")
            if paramfile.loopn == paramfile.loop:
                os.chdir('../Final/')
                scriptname = fm.writesubscript('abs',paramfile,clust)
                os.system(f"{clust.subco} {scriptname}")   
                with open(f'../{sumfilename}','a') as sumfile:
                    sumfile.write(f'Maximum number of loops reached\n')
            else:
                os.chdir('../')
                scriptname = fm.writesubscript('sparkplug',paramfile,clust)
                os.system(f"{clust.subco} {scriptname}")   
            with open(f'../{sumfilename}','a') as sumfile:
                sumfile.write(f'{emp:-^80}\n{emp:-^80}\n\n')
        else:
            with open(f'../{sumfilename}','a') as sumfile:
                sumfile.write("No new structure\n")
                sumfile.write(f'{emp:-^80}\n{emp:-^80}\n\n')
            os.chdir('../Final/')
            scriptname = fm.writesubscript('abs',paramfile,clust)
            os.system(f"{clust.subco} {scriptname}")   


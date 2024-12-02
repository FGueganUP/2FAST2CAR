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

# This script is used for final refining step
#
# In case of unconstrained CREST:
# It will first discard calculation that did not finish properly
# It will then submit frequency calculation for the most stable geometry
# It will then submit additional calculation cdft or NBO if selected
#
# In case of constrained CREST:
# It will then discard duplicate aka calculation with same first frequency and 
# same SCF energy at 6 digits
# It will then discard calculation with a first frequency not in 
# [- infinity, 0.7 first frequency of reference calculation]
# It will then submit additional calculation IRC

prevfilename = 'prevDFT.xyz'
loopfilename = 'loopDFT.xyz'

# Get cluster specifications
clust = rc.Cluster(fastcar_dir)
# Get keywords
kwds = rc.Keywords(fastcar_dir)
# Get information from parameters.txt file
paramfile = rp.Parameters('../parameters.tmp')
sumfilename = "summary.log"


now = datetime.datetime.now()
with open(f'../{sumfilename}','a') as sumfile:
    sumfile.write(f'Calling transmission at {now.time()}\n')

if paramfile.base_1 != paramfile.base_2:
    refine = True 
else:
    refine = False

if refine:
    newfile = []
# Analysis of TS
if paramfile.ref_freq:
    uniqueTS = []
    doublonTS = []
    garbageTS = []
    failedTS = []
    notcompletedTS = []

    
    if paramfile.check_st:
        loopfiles =  glob.glob(f'{paramfile.base_name_1}*_'
                               f'loop{paramfile.check_st}.out')
    else:
        loopfiles = glob.glob(f'{paramfile.base_name_1}*.out')

    # Sort TS by category
    #if paramfile.check_st:
    #    file_name = glob.glob('../Final/*_lowest.out')[0]
    #    startoutput = ro.Output(f'../Final/{file_name}')
    #else:
    #    startoutput = ro.Output(f'../{paramfile.experience_number}.out')
    #startoutput.readnormvec()
    #geor = np.zeros((startoutput.num_atom,3))
    #for i,at in enumerate(startoutput.coordinates.rstrip().split('\n')):
    #    for j,c in enumerate(at.split()[1:]):
    #        geor[i][j] = c
    goodfiles = list()
    for filename in loopfiles:
        outputfile = ro.Output(filename)
        if not outputfile.normterm:
            notcompletedTS.append((filename, 'Calculation failed'))
            if not 'failed' in filename:
                new_name = filename.split('.')[0] + '_failed.' + \
                           filename.split('.')[1]
                shutil.move(filename,new_name)
            continue
        outputfile.readfreqs() #Not converged?
        x = outputfile.freqs[0]
        if math.isnan(x):
            if not 'failed' in filename:
                new_name = filename.split('.')[0] + '_failed.' + \
                           filename.split('.')[1]
                shutil.move(filename,new_name)
            else:
                new_name = filename
            notcompletedTS.append((new_name, 'Freq calculation failed'))
            continue
        outputfile.readnrj()
        # Check if TS
        if x < 0:
            # Check if describe the good reaction
            outputfile.readnormvec() 
            #geoc = np.zeros((outputfile.num_atom,3))
            #for i,at in enumerate(outputfile.coordinates.rstrip()
            #                      .split('\n')):
            #    for j,c in enumerate(at.split()[1:]):
            #        geoc[i][j] = c
            #rot,rssd=sp.spatial.transform.Rotation.align_vectors(geor,geoc)
            #rot_mat = rot.as_matrix()
            #anormvec = list()
            #for v in outputfile.normvec:
            #    new_vec = rot_mat.dot(v[1])
            #    anormvec.append((v[0], new_vec))
            #isid = misc.compnv(startoutput.normvec,outputfile.normvec)
            #isid = misc.compnv(startoutput.normvec,anormvec)
            isid = misc.checknv(outputfile.normvec,paramfile.activ_ats)
            if isid:
                goodfiles.append(filename)
            else:
                if not 'garbage' in filename:
                    new_name = filename.split('.')[0] + '_garbage.' + \
                               filename.split('.')[1]
                    shutil.move(filename,new_name)
                else:
                    new_name = filename
                garbageTS.append((new_name, 'E SCF =',
                                  float(outputfile.last_SCF),
                                  'First frequency =', x))
        else:
            if not 'noTS' in filename:
                new_name = filename.split('.')[0] + '_noTS.' + \
                           filename.split('.')[1]
            else:
                new_name = filename
            shutil.move(filename,new_name)
            failedTS.append((new_name, 'E SCF =', float(outputfile.last_SCF), 
                            'First frequency =', x))
    rmsd_thresh = paramfile.RMSD_dft
    id_to_another = dict()
    newstrucs = list()
    doublons = list()
    if paramfile.check_st:
        prevfiles = glob.glob(f'{paramfile.base_name_1}*.out')
        prevfiles = [fn for fn in prevfiles if not f'loop{paramfile.check_st}'\
                     in fn and not 'doublon' in fn and not 'failed' in fn and \
                     not 'garbage' in fn and not 'noTS' in fn]
        fm.outstoxyzfile(prevfiles,prevfilename)
        prev_mols = io.loadallmols(prevfilename)
        fm.outstoxyzfile(goodfiles,loopfilename)
        loop_mols = io.loadallmols(loopfilename)
        geo_mols = io.loadallmols(loopfilename)
        sim = fm.geocomp(loop_mols,prev_mols,rmsd_thresh)
        nloop_mols = [cm for c,cm in enumerate(loop_mols) if c not in sim]
        ngeo_mols = [cm for c,cm in enumerate(geo_mols) if c not in sim]
        sn2 = [c for c,cm in enumerate(loop_mols) if c not in sim]
    else:
        fm.outstoxyzfile(goodfiles,loopfilename)
        nloop_mols = io.loadallmols(loopfilename)
        sn2 = range(0,len(nloop_mols))
    if len(nloop_mols) > 0:
        sim2 = fm.geocomp_onefile(nloop_mols,rmsd_thresh)
    else:
        sim2 = []
    if paramfile.check_st:
        for s in sim:
            filename = goodfiles[s]
            doublons.append(filename)
            if not 'doublon' in filename:
                outname = filename.split('.')[0] + '.out'
                new_name = outname.split('.')[0] + '_doublon.out'
                shutil.move(outname,new_name)
            else:
                new_name = filename
            id_to_another[new_name] = (prevfiles[sim[s][0]], sim[s][1])
    for s in sim2:
        filename = goodfiles[sn2[s]]
        doublons.append(filename)
        if not 'doublon' in filename:
            outname = filename.split('.')[0] + '.out'
            new_name = outname.split('.')[0] + '_doublon.out'
            shutil.move(outname,new_name)
        else:
            new_name = filename
        id_to_another[new_name] = (goodfiles[sn2[sim2[s][0]]], sim2[s][1])
    for fn in goodfiles:
        if fn not in doublons:
            newstrucs.append(fn)
        
    for filename in newstrucs:
        outputfile = ro.Output(filename)
        outputfile.readnrj()
        uniqueTS.append((filename, 'E SCF =',
                         float(outputfile.last_SCF),
                         'First frequency =', x))
        if refine:
            # Extract coordinate of 'unique' TS
            outputfile.extractcoord()
            # Create a new input file for large base calculation
            new_file_name = filename.replace(f'{base_name_1}', 
                                             f'{base_name_2}')
            new_file_name = new_file_name.replace('out', 'inp')
            soft = 'gaussian'
            kw = kwds.keywords[soft]['opt-ts']
            fm.writeg16inp(new_file_name,paramfile,12,
                           outputfile.coord_text,kw,refine=True)
            
    compname = f'{paramfile.base_name_1}*.out'
    # Check for identical structures
    #xyzlist = glob.glob(globname)
    #list_rmsd = list()
    #id_to_another = dict()
    #newstrucs = list()
    #sys.exit()
    #for filename in sorted(glob.glob(compname)):
    #    outputfile = ro.Output(filename)
    #    xyzname = filename.split('.')[0] + '.xyz'
    #    with open(xyzname,'w') as xyz:
    #        xyz.write(f'{outputfile.num_atom}\n\n{outputfile.coordinates}')
    #    if not paramfile.check_st or f'_loop{paramfile.check_st}' in filename:
    #        list_rmsd.append(xyzname)
    #for filename in list_rmsd:
    #    doublon,idfile,doub_rmsd = fm.xyzlistcleanerbis(filename, rmsd_thresh,
    #                                             compname.replace('out','xyz'))
    #    #doublon,idfile,doub_rmsd = fm.xyzlistcleanerbis(filename,0.5,
    #    #                                      compname.replace('out','xyz'))
    #    if doublon:
    #        id_to_another[filename.replace('xyz','out')] = (idfile.replace(
    #                                                    'xyz','out'),doub_rmsd)
    #    else:
    #        newstrucs.append(filename.replace('xyz','out'))
    if paramfile.check_st:
        pruname = f'pruningDFT_{paramfile.check_st}.log'
    else:
        pruname = f'pruningDFT.log'
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
            #sumfile.write(f'{len(id_to_another)} identical to another '
            #              'structure\n')
            heading = 'Structure identical to another one' 
            prunf.write(f'\n{heading:-^80}\n')
            for ita in id_to_another:
                prunf.write(f'{ita}: {id_to_another[ita]}\n')
    #for filename in sorted(glob.glob(compname.replace('out','xyz'))):
    #        os.remove(filename)
    for filename in id_to_another:
        #outputfile = ro.Output(filename.replace('xyz','out'))
        outputfile = ro.Output(filename)
        outputfile.readnrj()
        outputfile.readfreqs() 
        doublonTS.append(f'{filename} E SCF = {outputfile.last_SCF} First '
                         f'frequency = {outputfile.freqs[0]} id to ' 
                         f'{id_to_another[filename][0]} rmsd = '
                         f'{id_to_another[filename][1]}')
        if not 'doublon' in filename:
            outname = filename.split('.')[0] + '.out'
            new_name = outname.split('.')[0] + '_doublon.out'
            shutil.move(outname,new_name)

    #ADD REFINE CALCULATION SUBMITION
    uniqueTS.sort(key=lambda item: item[2])
    doublonTS.sort(key=lambda item: item[2])
    garbageTS.sort(key=lambda item: item[2])
    failedTS.sort(key=lambda item: item[2])
    if paramfile.check_st:
        wm = 'a'
    else:
        wm = 'w'
    fm.writelogTS(uniqueTS,doublonTS,garbageTS,failedTS,notcompletedTS,wm)
    with open(f'../{sumfilename}', 'a') as sumfile:
        sumfile.write('Informations about optimised TS written in '
                      'structures.log\n')
        sumfile.write(f'RMSD threshold: {paramfile.RMSD_dft}\n')
        sumfile.write(f'{len(uniqueTS)} new TS(s)\n') 
        sumfile.write(f'{len(doublonTS)} doublon(s)\n') 
        sumfile.write(f'{len(garbageTS)} TS(s) describing another reaction\n') 
        sumfile.write(f'{len(failedTS)} that are not TS\n') 
        sumfile.write(f'{len(notcompletedTS)} optimisation(s) failed\n') 


    if paramfile.check == 'none':
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
                sys.exit(0)
        os.chdir('../Final/')
        fm.writeabsscript(paramfile)
        os.system(f"sbatch script_abs.sub")   
    elif paramfile.check_st == 0:
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
                sys.exit(0)
        os.chdir('../')
        fm.writeloopscript(paramfile)
        os.system(f"sbatch script_loop.sub")   
    else:
        if len(uniqueTS) > 0:
            plname = glob.glob('../Final/*_lowest.out')[0]
            prevlowest = ro.Output(plname)
            prevlowest.readnrj()
            if uniqueTS[0][2] < prevlowest.last_SCF:
                lowestname = uniqueTS[0][0].split('.')[0] + '_lowest.' +\
                             uniqueTS[0][0].split('.')[1]
                shutil.copy(f'{uniqueTS[0][0]}', f'../Final/{lowestname}')
                os.remove(plname)
                with open(f'../{sumfilename}','a') as sumfile:
                    sumfile.write(f"New lowest energy structure: {lowestname}\n")
            else:
                with open(f'../{sumfilename}','a') as sumfile:
                    sumfile.write(f"No new lowest energy structure\n")
            if paramfile.check_st == paramfile.check:
                os.chdir('../Final/')
                fm.writeabsscript(paramfile)
                os.system(f"sbatch script_abs.sub")   
            else:
                os.chdir('../')
                fm.writeloopscript(paramfile)
                os.system(f"sbatch script_loop.sub")   
        else:
            with open(f'../{sumfilename}','a') as sumfile:
                sumfile.write("No new structure\n")
            os.chdir('../Final/')
            fm.writeabsscript(paramfile)
            os.system(f"sbatch script_abs.sub")   
   
# Analysis of geometry
else:
    E = []
    garbageGeo = []
    # Analysis of Geo to remove not converged calculation
    if paramfile.check_st:
        globname =  f'{paramfile.base_name_1}*_loop{paramfile.check_st}.out'
    else:
        globname =  f'{paramfile.base_name_1}*.out'
    for filename in sorted(glob.glob(globname)):
        outputfile = ro.Output(filename)
        outputfile.readnrj()
        if outputfile.normterm:
            E.append((filename, float(outputfile.last_SCF)))
        else:
            garbageGeo.append((filename, 'did not converged'))
    E.sort(key=lambda item: item[1])

    # Write in log file sorted geom by category
    if paramfile.check_st:
        wm = 'a'
    else:
        wm = 'w'
    fm.writelogmin(E,garbageGeo,wm)
    with open(f'../{sumfilename}', 'a') as sumfile:
        sumfile.write('Informations about optimised minima written in '
                      'structures.log\n')
        sumfile.write(f'{len(E)} new minima\n{len(garbageGeo)} not converged\n')

    if paramfile.check == 'none':
        if not paramfile.woc:
            if os.path.isdir('../Final/'):
                shutil.rmtree('../Final/')
            os.mkdir('../Final/')
        lowestname = E[0][0].split('.')[0] + '_lowest.' +\
                     E[0][0].split('.')[1]
        shutil.copy(E[0][0], f'../Final/{lowestname}')
        os.chdir('../Final/')
        fm.writeabsscript(paramfile)
        os.system(f"sbatch script_abs.sub")   
    elif paramfile.check_st == 0:
        if not paramfile.woc:
            if os.path.isdir('../Final/'):
                shutil.rmtree('../Final/')
            os.mkdir('../Final/')
        lowestname = E[0][0].split('.')[0] + '_lowest.' +\
                     E[0][0].split('.')[1]
        shutil.copy(E[0][0], f'../Final/{lowestname}')
        os.chdir('../')
        fm.writeloopscript(paramfile)
        os.system(f"sbatch script_loop.sub")   
    else:
        os.chdir('../')
        os.mkdir('COMP')
        nrjs = fm.readsummary('structures.log')
        outputfiles = dict()
        allmins = list()
        for lfn in nrjs:
            for ofn in nrjs[lfn]:
                if f'loop{paramfile.check_st}' in lfn:
                    ofn_comp = f'{ofn} check'
                else:
                    ofn_comp = f'{ofn} prev'
                outputfiles[ofn_comp] = ro.Output(f"Geo-CREST/{ofn}")
                outputfiles[ofn_comp].extractgeo()
                outputfiles[ofn_comp].readnrj()
        os.chdir('COMP/')
        for ofn in outputfiles:
            oxyzn_base = ofn.split('.')[0]
            oxyzn_ext = ofn.split()[1]
            oxyzn = f'{oxyzn_base}_{oxyzn_ext}.xyz'
            with open(oxyzn,'w') as of:
                of.write(f'{outputfiles[ofn].num_atom}\n')
                of.write(f'{outputfiles[ofn].last_SCF}\n')
                for item in outputfiles[ofn].coordinates:
                    of.write(' '.join(map(str, item)))
                    of.write('\n')

        corresp_check = dict()
        for ofn in outputfiles:
            oxyzn_base = ofn.split('.')[0]
            oxyzn_ext = ofn.split()[1]
            oxyzn = f'{oxyzn_base}_{oxyzn_ext}.xyz'
            with open(oxyzn,'r') as fout:
                flines = fout.readlines()
                nrj_recalc = float(flines[1])
            if not 'check' in ofn:
                allmins.append((ofn,nrj_recalc))
                continue
            result = subprocess.check_output(f'python -m spyrmsd -m {oxyzn} '
                                              '*prev.xyz', shell=True)
            prevfns = (sorted(glob.glob(f'*prev.xyz')))
            rmsd_values = result.split()
            rmsd_values = [float(item.decode()) for item in rmsd_values]
            mini = 999
            for r,rmsd in enumerate(rmsd_values):
                with open(prevfns[r],'r') as fout:
                    flines = fout.readlines()
                    nrj_prev = float(flines[1])
                nrjdiff = abs(nrj_prev-nrj_recalc)*627.51
                #if rmsd < paramfile.RMSD_dft_threshold and rmsd < mini and \
                #                                        nrjdiff < nrj_thresh:
                if rmsd < paramfile.RMSD_dft and rmsd < mini:
                    mini = rmsd
                    corresp_check[oxyzn]=(nrj_recalc,prevfns[r],rmsd,nrj_prev,
                                          nrjdiff)
                elif rmsd < mini:
                    mini = rmsd
                    closest = prevfns[r]
                    cl_nrj = nrj_prev
                    cl_diff = nrjdiff
                    cl_rmsd = rmsd
            if not oxyzn in corresp_check:
                corresp_check[oxyzn]=(f'new structure,{nrj_recalc},{closest},'
                                      f'{cl_rmsd},{cl_nrj},{cl_diff:.4f}')
                allmins.append((ofn,nrj_recalc))
        allmins.sort(key=lambda item: item[1])
        os.chdir('../')
        maxtitlesize = 0
        for elt in allmins:
            if len(elt[0]) > maxtitlesize:
                maxtitlesize = len(elt[0])
        maxtitlesize += 1
        with open('nrjs_check.log','w') as lf:
            for elt in allmins:
                lf.write(f'{elt[0]:{maxtitlesize}} {elt[1]}\n')
        with open('corresp_check.log','w') as lf:
            for oxyzn in corresp_check:
                lf.write(f'{oxyzn},{corresp_check[oxyzn]}\n')
        shutil.rmtree('COMP')
        if allmins[0][0].split()[1] == 'check':
            if not paramfile.woc:
                if os.path.isdir('Final/'):
                    shutil.rmtree('Final/')
                os.mkdir('Final/')
            lowestname = allmins[0][0].split()[0].split('.')[0] + '_lowest.' +\
                         allmins[0][0].split()[0].split('.')[1]
            shutil.copy(f'Geo-CREST/{allmins[0][0].split()[0]}', 
                        f'Final/{lowestname}')
        else:
            with open(sumfilename,'a') as sumfile:
                sumfile.write("No new lowest energy structure\n")
        if paramfile.check_st == paramfile.check:
            with open(sumfilename,'a') as sumfile:
                sumfile.write("Maximum number of loop reached\n")
            os.chdir('Final/')
            fm.writeabsscript(paramfile)
            os.system(f"sbatch script_abs.sub")   
        else:
            fm.writeloopscript(paramfile)
            os.system(f"sbatch script_loop.sub")   

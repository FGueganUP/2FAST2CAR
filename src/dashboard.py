#!/usr/bin/env python3
import glob
import os
import sys
import matplotlib                                                               
matplotlib.use('TkAgg')                                                         
import matplotlib.pyplot as plt                                                 
import numpy as np  
import re
import spyrmsd
from spyrmsd import io, rmsd
fastcar_dir = "/home/hgeindre/ProjectFASTCAR/FASTCAR2"
sys.path.append(fastcar_dir)
import filemanager as fm
import readoutput as ro
import readparam as rp
import misc
import plotfuncs as pf

options = ['Energies: CREST vs DFT','RMSD vs energies','CDFT','pruning',
           'Boltzmann']
conv = 627.509

ecs = fm.readepc()

paramfile = rp.Parameters('parameters.tmp')

if os.path.isdir('Geo-CREST'):
    calc_dir = 'Geo-CREST/'
elif os.path.isdir('TS-CREST'):
    calc_dir = 'TS-CREST/'
else:
    print('No directory with calculations')
    sys.exit(1)
print('Enter a number to chose an option:')
for o,opt in enumerate(options):
    print(f'-{o+1}: {opt}')
no = input() 
while True:
    try:
        no = int(no)-1
        if no < len(options) and no >= 0:
            break
        else:
            print('Please, enter an integer between 1 and '
                 f'{len(options)}')
            no = input() 
    except:
        print('Please, enter an integer')
        no = input() 

if no == 0:
    outfiles = dict()
    ls = list(filter(os.path.isdir, os.listdir()))
    crest_dirs = [d for d in ls if d.startswith('CREST')]
    crest_energies = dict()
    for crestd in crest_dirs:
        if '_' in crestd:
            ln = int(crestd.strip('CREST_')) 
        else:
            ln = 0
        crest_energies[ln] = list()
        with open(f'{crestd}/crest_conformers.xyz','r' ) as confile:
            conflines = confile.readlines()
            nat = int(conflines[0])
            cn = len(conflines) // (nat+2) 
            l = 1
            for c in range(0,cn):
                nrj = float(conflines[l])
                l += nat + 2
                crest_energies[ln].append(nrj)
                 
    xs = dict()
    ys = dict()
    xs['all'] = list()
    ys['all'] = list()
    lowestcrest = 0
    lowestdft = 0
    for filename in sorted(glob.glob(f'{calc_dir}'
                                     f'*{paramfile.experience_number}*.out')):
        if 'doublon' in filename or 'garbage' in filename \
                                  or 'failed' in filename:
            continue
        if 'loop' in filename:
            ln = int(filename.split('.')[0].split('_')[-1].strip('loop'))
            cn = int(filename.split('.')[0].split('-')[-1].split('_')[0])
        else:
            ln = 0
            cn = int(filename.split('.')[0].split('-')[-1])
        outfiles[filename] = ro.Output(filename)
        outfiles[filename].readnrj()
        if not ln in xs:
            xs[ln] = list()
            ys[ln] = list()
        if crest_energies[ln][cn-1] < lowestcrest:
            lowestcrest = crest_energies[ln][cn-1]
        if outfiles[filename].last_SCF < lowestdft:
            lowestdft = outfiles[filename].last_SCF
        xs[ln].append(crest_energies[ln][cn-1])
        ys[ln].append(outfiles[filename].last_SCF)
        xs['all'].append(crest_energies[ln][cn-1])
        ys['all'].append(outfiles[filename].last_SCF)
        
    for ln in xs:
        xs[ln] = [(x-lowestcrest)*conv for x in xs[ln]]
        ys[ln] = [(y-lowestdft)*conv for y in ys[ln]]
    r = np.corrcoef(xs['all'],ys['all'])                                                          
    
    lgd = list()
    plt.title("Energy correlation")                                                 
    plt.xlabel("CREST energies (kcal/mol)")
    plt.ylabel("DFT energies (kcal/mol)")
    colorlist = list(matplotlib.colors.TABLEAU_COLORS.keys())
    for ln in xs:
        if ln == 'all':
            continue
        color = colorlist[ln]
        plt.scatter(xs[ln],ys[ln],c=matplotlib.colors.TABLEAU_COLORS[color])
        lgd.append(f'loop {ln}')
    plt.legend(lgd)
    plt.annotate(f'R = {r[0,1]}', xy=(0.05, 0.95), xycoords='axes fraction')
    plt.show()                                                                      
    plt.close() 
elif no == 1:
    data = dict()
    for filename in sorted(glob.glob(f'{calc_dir}'
                                     f'*{paramfile.experience_number}*.out')):
        if 'doublon' in filename or 'garbage' in filename \
                                  or 'failed' in filename:
        #if 'garbage' in filename or 'failed' in filename:
            continue
        outputfile = ro.Output(filename)
        if not outputfile.normterm:
            continue
        outputfile.readnrj()
        fm.outstoxyzfile([filename],'tmp.xyz')
        mol = io.loadmol('tmp.xyz')
        mol.strip()
        data[filename] = (outputfile.last_SCF,mol)
    xs = list()
    ys = list()
    for f,filename in enumerate(data):
        coords1 = data[filename][1].coordinates
        anum1 = data[filename][1].atomicnums
        adj1 = data[filename][1].adjacency_matrix
        for f2,filename2 in enumerate(data):
            if f2 >= f:
                continue
            coords2 = data[filename2][1].coordinates
            anum2 = data[filename2][1].atomicnums
            adj2 = data[filename2][1].adjacency_matrix
            rmsd_value = rmsd.symmrmsd(coords1,coords2,anum1,anum2,adj1,adj2, 
                                       minimize=True)
            deltaE = abs(data[filename][0] - data[filename2][0])*conv
            xs.append(deltaE)
            ys.append(rmsd_value)
            if rmsd_value < 0.3:
                print(filename,filename2,rmsd_value)
    plt.title("Energy difference and rmsd correlation")
    plt.xlabel("Energy difference (kcal/mol)")
    plt.ylabel("RMSD")
    plt.xscale("log")
    plt.scatter(xs,ys)
    plt.show()                                                                      
    plt.close() 
    
elif no == 2:
    cdftfiles = list()
    for filename in sorted(glob.glob(f'Final/*{paramfile.experience_number}*'
                                      '.out')):
        outputfile = ro.Output(filename)
        outputfile.readnrj()
        # Extract NPA
        if 'cdft' in filename:
            outputfile.readcharges()
            cdftfiles.append(outputfile)
    # Calculate some conceptual dft descriptors
    if len(cdftfiles) > 0:
        cdft_gds = misc.calc_cdft_gds(cdftfiles)
        cdft_funcs = misc.calc_cdft_funcs(cdftfiles,cdft_gds)
        fm.writecdftlog(cdft_gds,cdft_funcs,cdftfiles)
        pf.plot_cdft(cdftfiles[0],cdft_funcs,ecs)

elif no == 3:
    pruned = dict()
    for filename in sorted(glob.glob(f'{calc_dir}pruningCREST*.log')):
        with open(filename,'r') as prunfile:
            prunlines = prunfile.readlines()
            read = False
            for line in prunlines:
                if not read and 'Structure identical' in line:
                    read = True
                elif read and '----------' in line or line.isspace():
                    read = False
                elif read:
                    sline = line.split()
                    pruned[sline[0].strip(':')] = (sline[1].strip('(\','),
                                             float(sline[2].strip(')\'')))
    print('Which structure to show?')
    spruned = {k: v for k, v in sorted(pruned.items(), 
               key=lambda item: item[1][1], reverse = True)}
    for r,rej in enumerate(spruned):
        print(f'{r+1}: {rej} {spruned[rej]}')
    nr = int(input())
    rstruc = list(spruned.keys())[nr-1]
    istruc = spruned[rstruc][0]
    rmsd = spruned[rstruc][1]
    rmatch = re.findall(r'\d+', rstruc)
    imatch = re.findall(r'\d+', istruc)
    rnr = int(rmatch[-1])
    rni = int(imatch[-1])
    if len(rmatch) > 1:
        rloop = int(rmatch[0])
        rconfile = f'{calc_dir}crest_conformers_{rloop}.xyz'
    else:
        rconfile = f'{calc_dir}crest_conformers.xyz'
    if len(imatch) > 1:
        iloop = int(imatch[0])
        ionfile = f'{calc_dir}crest_conformers_{iloop}.xyz'
    else:
        iconfile = f'{calc_dir}crest_conformers.xyz'
    rejgeo = fm.xyzfinder(rconfile,rnr-1)
    idgeo = fm.xyzfinder(iconfile,rni-1)
    pf.superposeGeos(rejgeo,idgeo,rmsd,ecs)

elif no == 4:
    ewin = 3
    ls = list(filter(os.path.isdir, os.listdir()))
    crest_dirs = [d for d in ls if d.startswith('CREST')]
    crest_energies = dict()
    min_loop = dict()
    for crestd in crest_dirs:
        if '_' in crestd:
            ln = int(crestd.strip('CREST_')) 
        else:
            ln = 0
        if not ln in min_loop:
            min_loop[ln] = 0
        crest_energies[ln] = list()
        with open(f'{crestd}/crest_conformers.xyz','r' ) as confile:
            conflines = confile.readlines()
            nat = int(conflines[0])
            cn = len(conflines) // (nat+2) 
            l = 1
            for c in range(0,cn):
                nrj = float(conflines[l])
                l += nat + 2
                crest_energies[ln].append(nrj)
                if nrj < min_loop[ln]:
                   min_loop[ln] = nrj 
                

    mini = 0
    nrjs = dict()
    crest_nrj = dict()
    for filename in sorted(glob.glob(f'{calc_dir}'
                                     f'*{paramfile.experience_number}*.out')):
        if 'doublon' in filename or 'garbage' in filename \
                                  or 'failed' in filename:
            continue
        outputfile = ro.Output(filename)
        if not outputfile.normterm:
            continue
        outputfile.readnrj()
        nrjs[filename] = outputfile.last_SCF
        if outputfile.last_SCF < mini:
            mini = outputfile.last_SCF
            lowestfile = filename
    exp = dict()
    exp_sum = 0
    pops_loop = dict()
    for filename in nrjs:
        exp[filename] = np.exp(-(nrjs[filename]-nrjs[lowestfile])\
                        /(273*0.0000031668))
        exp_sum += exp[filename]
    for filename in exp:
        if 'loop' in filename:
            ln = int(filename.split('.')[0].split('_')[-1].strip('loop'))
            cn = int(filename.split('.')[0].split('-')[-1].split('_')[0])
        else:
            ln = 0
            cn = int(filename.split('.')[0].split('-')[-1])
        deltaE_crest = (crest_energies[ln][cn]-min_loop[ln])*conv
        if deltaE_crest > ewin:
            continue
        pop = (exp[filename]/exp_sum)*100
        if not ln in pops_loop:
            pops_loop[ln] = 0
        pops_loop[ln] += pop
    x = list()
    lgd = list()
    pop_tot = 0
    for ln in sorted(list(pops_loop.keys())):
        x.append(pops_loop[ln])
        lgd.append(f'loop {ln}: {pops_loop[ln]:1.1f}%')
        pop_tot += pops_loop[ln]
    #plt.pie(x,labels=lgd, autopct='%1.1f%%', startangle=90)
    plt.pie(x,startangle=90)
    #plt.hist(x,label=lgd)
    plt.legend(lgd)
    plt.tight_layout()
    plt.show()


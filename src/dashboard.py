#!/usr/bin/env python3
import glob
import os
import sys
import datetime
from tqdm import tqdm
import matplotlib                                                               
matplotlib.use('TkAgg')                                                         
import matplotlib.pyplot as plt                                                 
import numpy as np  
import re
import spyrmsd
from spyrmsd import io, rmsd
from textwrap import wrap
fastcar_dir = "/home/hgeindre/ProjectFASTCAR/FASTCAR2"
sys.path.append(fastcar_dir)
import filemanager as fm
import readoutput as ro
import readparam as rp
import misc
import plotfuncs as pf

#options = ['OPT energies vs CREST energies','OPT RMSD vs OPT energies',
#           'OPT RMSD vs CREST RMSD','CDFT','Pruning','Boltzmann','Structures',
#           'Write energies csv','Execute command']
options = ['OPT energies vs CREST energies','OPT RMSD vs OPT energies',
           'OPT RMSD vs CREST RMSD','CDFT','Pruning','Boltzmann','Structures',
           'Write energies csv','Time stats']
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
    xs[-1] = list()
    ys[-1] = list()
    lowestcrest = 0
    lowestopt = 0
    start_nrj = 'none'
    for filename in sorted(glob.glob(f'{calc_dir}'
                                     f'*{paramfile.calc_name}*.out')):
        if 'duplicate' in filename or 'garbage' in filename \
                                 or 'failed' in filename or 'noTS' in filename:
            continue
        if 'loop' in filename:
            ln = int(filename.split('.')[0].split('_')[-1].strip('loop'))
            cn = int(filename.split('.')[0].split('-')[-1].split('_')[0])
        else:
            ln = 0
            cn = int(filename.split('.')[0].split('-')[-1])
        outfiles[filename] = ro.Output(filename)
        if not outfiles[filename].normterm:
            continue
        if cn == 0:
            start_nrj = outfiles[filename].last_SCF
            if start_nrj < lowestopt:
                lowestopt = start_nrj
            continue
        if not ln in xs:
            xs[ln] = list()
            ys[ln] = list()
        if crest_energies[ln][cn-1] < lowestcrest:
            lowestcrest = crest_energies[ln][cn-1]
        if outfiles[filename].last_SCF < lowestopt:
            lowestopt = outfiles[filename].last_SCF
        xs[ln].append(crest_energies[ln][cn-1])
        ys[ln].append(outfiles[filename].last_SCF)
        xs[-1].append(crest_energies[ln][cn-1])
        ys[-1].append(outfiles[filename].last_SCF)
        
    for ln in xs:
        xs[ln] = [(x-lowestcrest)*conv for x in xs[ln]]
        ys[ln] = [(y-lowestopt)*conv for y in ys[ln]]
    if start_nrj != 'none':
        start_nrj = (start_nrj-lowestopt)*conv
    r = np.corrcoef(xs[-1],ys[-1])                                                          
    
    lgd = list()
    plt.rcParams["figure.autolayout"] = True
    plt.rcParams["figure.figsize"] = [15, 10]
    matplotlib.rcParams.update({'font.size': 22})
    plt.title("Energy correlation")                                                 
    plt.xlabel("CREST energies (kcal/mol)")
    plt.ylabel("OPT energies (kcal/mol)")
    colorlist = list(matplotlib.colors.TABLEAU_COLORS.keys())
    for ln in sorted(xs):
        if ln == -1:
            continue
        if ln+1 > len(colorlist):
            lnc = (ln+1) % len(colorlist) - 1
        else:
            lnc = ln
        color = colorlist[lnc]
        plt.scatter(xs[ln],ys[ln],c=matplotlib.colors.TABLEAU_COLORS[color])
        lgd.append(f'loop {ln}')
    if start_nrj != 'none':
        plt.axhline(y=start_nrj, color='r', linestyle='-')
    plt.legend(lgd)
    plt.annotate(f'R = {r[0,1]}', xy=(0.25, 0.95), xycoords='axes fraction')
    plt.show()                                                                      
    plt.close() 
elif no == 1:
    data = dict()
    for filename in sorted(glob.glob(f'{calc_dir}'
                                     f'*{paramfile.calc_name}*.out')):
        if 'duplicate' in filename or 'garbage' in filename \
                                 or 'failed' in filename or 'noTS' in filename:
        #if 'garbage' in filename or 'failed' in filename:
            continue
        outputfile = ro.Output(filename)
        if not outputfile.normterm:
            continue
        fm.outstoxyzfile([filename],'tmp.xyz')
        mol = io.loadmol('tmp.xyz')
        mol.strip()
        data[filename] = (outputfile.last_SCF,mol)
    os.remove('tmp.xyz')
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
    plt.rcParams["figure.autolayout"] = True
    plt.rcParams["figure.figsize"] = [15, 10]
    matplotlib.rcParams.update({'font.size': 22})
    plt.title("Energy difference and rmsd correlation")
    plt.xlabel("Energy difference (kcal/mol)")
    plt.ylabel("RMSD")
    plt.xscale("log")
    plt.scatter(xs,ys)
    plt.show()                                                                      
    plt.close() 
    
elif no == 2:
    data = dict()
    for filename in sorted(glob.glob(f'{calc_dir}'
                                     f'*{paramfile.calc_name}*.out')):
        if 'garbage' in filename or 'failed' in filename or 'noTS' in filename:
            continue
        outputfile = ro.Output(filename)
        if not outputfile.normterm:
            continue
        if 'loop' in filename:
            if 'duplicate' in filename:
                ln = int(filename.split('.')[0].split('_')[-2].strip('loop'))
            else:
                ln = int(filename.split('.')[0].split('_')[-1].strip('loop'))
            cn = int(filename.split('.')[0].split('-')[-1].split('_')[0])
            ccfname = f'crest_conformers_{ln}.xyz'
        else:
            ln = 0
            cn = int(filename.split('.')[0].split('-')[-1].split('_')[0])
            if cn == 0:
                continue
            ccfname = f'crest_conformers.xyz'
        fm.outstoxyzfile([filename],'tmp.xyz')
        opt_mol = io.loadmol('tmp.xyz')
        opt_mol.strip()
        crest_geo, crest_atnums = fm.xyzfinder(f'{calc_dir}{ccfname}',cn)
        fm.writexyzfile(crest_geo,crest_atnums,'tmp.xyz')
        crest_mol = io.loadmol('tmp.xyz')
        crest_mol.strip()
        data[filename] = (crest_mol,opt_mol)
    os.remove('tmp.xyz')
    xs = list()
    ys = list()
    for f,filename in enumerate(data):
        crest_coords1 = data[filename][0].coordinates
        crest_anum1 = data[filename][0].atomicnums
        crest_adj1 = data[filename][0].adjacency_matrix
        opt_coords1 = data[filename][1].coordinates
        opt_anum1 = data[filename][1].atomicnums
        opt_adj1 = data[filename][1].adjacency_matrix
        for f2,filename2 in enumerate(data):
            if f2 >= f:
                continue
            crest_coords2 = data[filename2][0].coordinates
            crest_anum2 = data[filename2][0].atomicnums
            crest_adj2 = data[filename2][0].adjacency_matrix
            crest_rmsd_value = rmsd.symmrmsd(crest_coords1,crest_coords2,
                                             crest_anum1,crest_anum2,
                                             crest_adj1,crest_adj2,
                                             minimize=True)
            opt_coords2 = data[filename2][1].coordinates
            opt_anum2 = data[filename2][1].atomicnums
            opt_adj2 = data[filename2][1].adjacency_matrix
            opt_rmsd_value = rmsd.symmrmsd(opt_coords1,opt_coords2,
                                           opt_anum1,opt_anum2,
                                           opt_adj1,opt_adj2,minimize=True)
            xs.append(crest_rmsd_value)
            ys.append(opt_rmsd_value)
    plt.rcParams["figure.autolayout"] = True
    plt.rcParams["figure.figsize"] = [15, 10]
    matplotlib.rcParams.update({'font.size': 22})
    plt.title("RMSD correlation")
    plt.xlabel("RMSD CREST")
    plt.ylabel("RMSD OPT")
    plt.scatter(xs,ys)
    plt.show()                                                                      
    plt.close() 
    
elif no == 3:
    cdftfiles = list()
    for filename in sorted(glob.glob(f'Final/*{paramfile.calc_name}*'
                                      '.out')):
        outputfile = ro.Output(filename)
        # Extract NPA
        if 'cdft' in filename:
            cdftfiles.append(outputfile)
    # Calculate some conceptual dft descriptors
    if len(cdftfiles) > 0:
        cdft_gds = misc.calc_cdft_gds(cdftfiles)
        cdft_funcs = misc.calc_cdft_funcs(cdftfiles,cdft_gds)
        fm.writecdftlog(cdft_gds,cdft_funcs,cdftfiles)
        pf.plot_cdft(cdftfiles[0],cdft_funcs,ecs)

elif no == 4:
    pruned = dict()
    for filename in sorted(glob.glob(f'{calc_dir}pruningCREST*.log')):
        if '_' in filename:
            ln=int(filename.replace(calc_dir,'').split('.')[0].split('_')[-1])
        else:
            ln=0
        with open(filename,'r') as prunfile:
            prunlines = prunfile.readlines()
            read = False
            for line in prunlines:
                if not read and 'Structure identical' in line:
                    read = True
                    if 'previous' in line:
                        prev = True
                    else:
                        prev = False
                elif read and '----------' in line or line.isspace():
                    read = False
                elif read:
                    sline = line.split()
                    cn = sline[0].strip(':')
                    if prev:
                        if 'loop' in sline[1]:
                            ln2 = int(sline[1].split('.')[0].split('_')[-1]
                                              .strip('loop'))
                            cn2 = int(sline[1].split('.')[0].split('-')[-1]
                                              .split('_')[0])
                        else:
                            ln2 = 0
                            cn2 = int(sline[1].split('.')[0].split('-')[-1])
                    else:
                        ln2 = ln
                        cn2 =  int(sline[1])
                    pruned[f'{cn} {ln}'] = (f'{cn2} {ln2}',
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
    rnr = int(rstruc.split()[0])
    rloop = int(rstruc.split()[1])
    rni = int(istruc.split()[0])
    iloop = int(istruc.split()[1])
    if rloop:
        rconfile = f'{calc_dir}crest_conformers_{rloop}.xyz'
    else:
        rconfile = f'{calc_dir}crest_conformers.xyz'
    if iloop:
        iconfile = f'{calc_dir}crest_conformers_{iloop}.xyz'
    else:
        iconfile = f'{calc_dir}crest_conformers.xyz'
    rejgeo, rejatnums = fm.xyzfinder(rconfile,rnr)
    idgeo, idatnums = fm.xyzfinder(iconfile,rni)
    pf.superposeGeos(rejgeo,rejatnums,idgeo,idatnums,rmsd,ecs)

elif no == 5:
    #ewin = 3
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
                                     f'*{paramfile.calc_name}*.out')):
        if 'duplicate' in filename or 'garbage' in filename \
                                 or 'failed' in filename or 'noTS' in filename:
            continue
        outputfile = ro.Output(filename)
        if not outputfile.normterm:
            continue
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
        #if deltaE_crest > ewin:
        #    continue
        pop = (exp[filename]/exp_sum)*100
        if not ln in pops_loop:
            pops_loop[ln] = 0
        pops_loop[ln] += pop
    x = list()
    lgd = list()
    pop_tot = 0
    plt.rcParams["figure.autolayout"] = True
    plt.rcParams["figure.figsize"] = [15, 10]
    matplotlib.rcParams.update({'font.size': 18})
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

elif no == 6:
    # Rajouter minima
    steps = list()
    counts = dict()
    #counts['CREST conformers'] = list()
    counts['CREST new conformers'] = list()
    counts['CREST duplicates'] = list()
    counts['CREST duplicates previous loop'] = list()
    counts['OPT new TS'] = list()
    counts['OPT duplicates'] = list()
    counts['OPT garbage'] = list()
    counts['OPT no TS'] = list()
    counts['OPT failed'] = list()
    with open('summary.log', 'r') as sumfile:
        sumlines = sumfile.readlines()
        for line in sumlines:
            if '-------------Loop' in line:
                ln = int(line.replace('-','').split()[1])
                if ln > 0:
                    #counts['CREST conformers'].append(crest_tc)
                    #counts['CREST conformers'].append(0)
                    counts['CREST new conformers'].append(crest_nc)
                    counts['CREST new conformers'].append(0)
                    counts['CREST duplicates'].append(crest_d)
                    counts['CREST duplicates'].append(0)
                    counts['CREST duplicates previous loop'].append(crest_pd)
                    counts['CREST duplicates previous loop'].append(0)
                    counts['OPT new TS'].append(0)
                    counts['OPT new TS'].append(opt_nts)
                    counts['OPT duplicates'].append(0)
                    counts['OPT duplicates'].append(opt_d)
                    counts['OPT garbage'].append(0)
                    counts['OPT garbage'].append(opt_g)
                    counts['OPT no TS'].append(0)
                    counts['OPT no TS'].append(opt_nots)
                    counts['OPT failed'].append(0)
                    counts['OPT failed'].append(opt_f)
                steps.append(f'loop {ln}\n CREST')
                steps.append(f'loop {ln}\n OPT')
                crest_tc = 0
                crest_nc = 0
                crest_d = 0
                crest_pd = 0
                opt_nts = 0
                opt_d = 0
                opt_g = 0
                opt_nots = 0
                opt_f = 0
            elif 'conformers found by CREST' in line:
                crest_tc = int(line.split()[0])
            elif 'new conformers' in line:
                crest_nc = int(line.split()[0])
            elif 'identical to another structure from the same loop' in line:
                crest_d = int(line.split()[0])
            elif 'identical to a structure from a previous loop' in line:
                crest_pd = int(line.split()[0])
            elif 'new TS(s)' in line:
                opt_nts = int(line.split()[0])
            elif 'duplicate(s)' in line:
                opt_d = int(line.split()[0])
            elif 'TS(s) describing another reaction' in line:
                opt_g = int(line.split()[0])
            elif 'that are not TS' in line:
                opt_nots = int(line.split()[0])
            elif 'optimisation(s) failed' in line:
                opt_f = int(line.split()[0])
        #counts['CREST conformers'].append(crest_tc)
        #counts['CREST conformers'].append(0)
        counts['CREST new conformers'].append(crest_nc)
        counts['CREST new conformers'].append(0)
        counts['CREST duplicates'].append(crest_d)
        counts['CREST duplicates'].append(0)
        counts['CREST duplicates previous loop'].append(crest_pd)
        counts['CREST duplicates previous loop'].append(0)
        counts['OPT new TS'].append(0)
        counts['OPT new TS'].append(opt_nts)
        counts['OPT duplicates'].append(0)
        counts['OPT duplicates'].append(opt_d)
        counts['OPT garbage'].append(0)
        counts['OPT garbage'].append(opt_g)
        counts['OPT no TS'].append(0)
        counts['OPT no TS'].append(opt_nots)
        counts['OPT failed'].append(0)
        counts['OPT failed'].append(opt_f)

    plt.rcParams["figure.autolayout"] = True
    plt.rcParams["figure.figsize"] = [15, 10]
    matplotlib.rcParams.update({'font.size': 18})
    fig, ax = plt.subplots()
    bottom = np.zeros(len(steps))
    width = 0.5
    
    for boolean, number in counts.items():
        p = ax.bar(steps, number, width, label='\n'.join(wrap(boolean, 15)), 
                   bottom=bottom)
        bottom += number
    
    #ax.set_title("Number of penguins with above average body mass")
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    for c in ax.containers:
        #labels = [v.get_height() if v.get_height() > 0 else '' for v in c]
        #ax.bar_label(c, labels=labels, label_type='center')
        ax.bar_label(c, fmt=lambda x: f'{x:.0f}' if x > 0 else '', 
                     label_type='center')
    
    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    #ax.legend(loc="upper right")
    
    plt.show()
elif no == 7:
    lowestname = glob.glob('Final/*_lowest.out')
    file_name = lowestname[0]
    lowestfile = ro.Output(file_name) 
    lowestnrj = lowestfile.last_SCF
    data = list()
    for filename in sorted(glob.glob(f'{calc_dir}'
                                     f'*{paramfile.calc_name}*.out')):
        if 'duplicate' in filename or 'garbage' in filename or 'failed' in \
                          filename or 'noTS' in filename:
            continue
        outputfile = ro.Output(filename)
        if not outputfile.normterm:
            continue
        e = outputfile.last_SCF
        deltaE = (outputfile.last_SCF - lowestnrj)*conv
        data.append((filename.split('/')[-1],e,deltaE))
    data.sort(key=lambda k: k[2]) 
    csvname = f'{paramfile.calc_name}_nrjs.csv'
    with open(csvname,'w') as csv:
        csv.write(f'{paramfile.calc_name}\n')
        csv.write(f'Name;E(ua);deltaE(kcal/mol)\n')
        for d in data:
            csv.write(f'{d[0]};{d[1]};{d[2]}\n')
    print(f'File {csvname} created')
elif no == 8:
    trt = datetime.timedelta()
    byloop = dict()
    bytype = dict()
    bytype['success'] = datetime.timedelta()
    bytype['duplicate'] = datetime.timedelta()
    bytype['garbage'] = datetime.timedelta()
    bytype['failed'] = datetime.timedelta()
    bytype['noTS'] = datetime.timedelta()
    nopt = 0 #Number of optimisations
    nwt = 0 #Number of failed calculations without time
    with tqdm(total=len(glob.glob(f'{calc_dir}'f'*{paramfile.calc_name}*.out')
                       )) as bar:
        for filename in sorted(glob.glob(f'{calc_dir}'
                                         f'*{paramfile.calc_name}*.out')):
            if 'loop' in filename:
                loop = re.search(r'loop(\d+)',filename)
                nloop = int(loop.group(1))
            else:
                nloop = 0
            outputfile = ro.Output(filename) 
            if outputfile.cputime != 'none':
                for rt in outputfile.cputime:
                    trt += rt
                nopt += 1
            else:
                nwt += 1
            if 'duplicate' in filename or 'doublon' in filename:
                for rt in outputfile.cputime:
                    bytype['duplicate'] += rt
            elif 'garbage' in filename:
                for rt in outputfile.cputime:
                    bytype['garbage'] += rt
            elif 'failed' in filename:
                if outputfile.cputime != 'none':
                    for rt in outputfile.cputime:
                        bytype['failed'] += rt
            elif 'noTS' in filename:
                for rt in outputfile.cputime:
                    bytype['noTS'] += rt
            else:
                if outputfile.cputime != 'none':
                    for rt in outputfile.cputime:
                        bytype['success'] += rt
            if not nloop in byloop:
                byloop[nloop] = datetime.timedelta()
            for rt in outputfile.cputime:
                if outputfile.cputime != 'none':
                    byloop[nloop] += rt
            bar.update(1)
    av = trt / nopt
    print(trt,nopt,av,nwt)
    print(byloop)
    print(bytype)
    print(nopt,trt/datetime.timedelta(hours=1))
    
    

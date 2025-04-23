#!/usr/bin/env python3
import os
import sys
import glob
import datetime
fastcar_dir = "/home/hgeindre/ProjectFASTCAR/FASTCAR2"
sys.path.append(fastcar_dir)
import readparam as rp
import readoutput as ro
import filemanager as fm
import readconfig as rc
import misc

sumfilename = "summary.log"
now = datetime.datetime.now()
date = now.strftime("%d/%m/%y %H:%M")
with open(f"../{sumfilename}",'a') as sumfile:
    sumfile.write(f'Calling abs.py {date}\n')


# Get cluster specifications
clust = rc.Cluster(fastcar_dir)
# Get keywords
kwds = rc.Keywords(fastcar_dir)
# Get information from parameters.txt file
paramfile = rp.Parameters('../parameters.tmp')
paramfile.getopt()
paramfile.getwoc()
paramfile.getadcalc()

if paramfile.adcalc == 'none':
    # Call wheels.py
    if paramfile.woc:
        dep = False
    else:
        dep = True
    scriptname = fm.writesubscript('wheels',paramfile,clust,dependency=dep)
    os.system(f'{clust.subco} {scriptname}')
    sys.exit(0)

filenames = dict()
file_name = dict()
if 'all' in paramfile.scope:
    filenames['all'] = list()
    if os.path.isdir('../Geo-CREST'):
        calc_dir = 'Geo-CREST/'
    elif os.path.isdir('../TS-CREST'):
        calc_dir = 'TS-CREST/'
    else:
        print('No directory with calculations')
        sys.exit(1)
    for filename in sorted(glob.glob(f'../{calc_dir}'
                                     f'*{paramfile.calc_name}*.out')):
        if 'duplicate' in filename or 'garbage' in filename \
                                 or 'failed' in filename or 'noTS' in filename:
            continue
        filenames['all'].append(filename)
    file_name['all'] = "all conformers"
elif 'lowest' in paramfile.scope:
    filenames['lowest'] = list()
    filenames['lowest'].append(glob.glob('*_lowest.out')[0])
    file_name['lowest'] = filenames['lowest'][0]
else:
    print(f'The value of the scope variable is abnormal: {paramfile.scope}')
    sys.exit(1)

for a,ac in enumerate(paramfile.adcalc):
    # Launch IRC calculation with the most suitable TS if chosen
    #soft = 'gaussian'
    if a >= len(paramfile.ad_software):
        software = paramfile.ad_software[-1]
    else:
        software = paramfile.ad_software[a]
    kw = (kwds.keywords[software][ac],ac)
    if a >= len(paramfile.scope):
        scope = paramfile.scope[-1]
    else:
        scope = paramfile.scope[a]
    if a >= len(paramfile.ad_method):
        method = paramfile.ad_method[-1]
    else:
        method = paramfile.ad_method[a]
    if a >= len(paramfile.ad_basis_set):
        basis = paramfile.ad_basis_set[-1]
    else:
        basis = paramfile.ad_basis_set[a]
    if a >= len(paramfile.ad_solvent):
        solvent = paramfile.ad_solvent[-1]
    else:
        solvent = paramfile.ad_solvent[a]
    if a >= len(paramfile.ad_dispersion):
        disp = paramfile.ad_dispersion[-1]
    else:
        disp = paramfile.ad_dispersion[a]
    gfn = file_name[scope]
    if ac == 'single-point' or ac == 'wfx':
        for fn in filenames[scope]:
            outfile = ro.Output(fn) 
            jobname = f'{ac}_{outfile.base_name}'
            inpname = f'{jobname}.inp'
            fm.writeinp(inpname,kw,outfile.coordinates,outfile.atnums,
                        paramfile,method,basis,paramfile.charge,
                        paramfile.mult,paramfile.cpus,solvent,disp,
                        soft=software)
            fm.writecalcsubscript(paramfile,jobname,clust,soft=software)
            if not paramfile.woc:
                os.system(f"{clust.subco} {jobname}.sub")
        with open(f'../{sumfilename}','a') as sumfile:
            sumfile.write(f'Single point calculations submitted for {gfn}\n')
    if ac == 'reopt':
        for fn in filenames[scope]:
            outfile = ro.Output(fn) 
            jobname = f'{ac}_{outfile.base_name}'
            inpname = f'{jobname}.inp'
            #fm.writeg16inp(inpname,paramfile,12,outfile.coordinates,
            #               outfile.atnums,kw,ac=a)
            fm.writeinp(inpname,kw,outfile.coordinates,outfile.atnums,
                        paramfile,method,basis,paramfile.charge,
                        paramfile.mult,paramfile.cpus,solvent,disp,
                        soft=software)
            fm.writecalcsubscript(paramfile,jobname,clust,soft=software)
            if not paramfile.woc:
                os.system(f"{clust.subco} {jobname}.sub")
        with open(f'../{sumfilename}','a') as sumfile:
            sumfile.write(f'Reoptimisation calculations submitted for {gfn}\n')
    elif ac == 'irc':
        for fn in filenames[scope]:
            outfile = ro.Output(fn) 
            IRC = ["forward", "reverse"]
            for item in IRC:
                jobname = f'{str(item)}_{ac}_{outfile.base_name}'
                inpname = f'{jobname}.inp'
                chkname = f'{jobname}.chk'
                kw = (kw[0].replace('direction',item),kw[1])
                kw2 = kwds.keywords[software]['irc-opt']
                fm.writeinp(inpname,kw,outfile.coordinates,outfile.atnums,
                            paramfile,method,basis,paramfile.charge,
                            paramfile.mult,paramfile.cpus,solvent,disp,
                            soft=software,chk=chkname,supcalc=True,
                            keywords2=kw2)
                #fm.writeg16inp(inpname,paramfile,12,outfile.coordinates,
                #              outfile.atnums,kw,chkn=chkname,link1=True,
                #               keywords2=kw2,ac=a)
                fm.writecalcsubscript(paramfile,jobname,clust,
                                      time='10-00:00:00',soft=software)
                if not paramfile.woc:
                    os.system(f"{clust.subco} {jobname}.sub")
        with open(f'../{sumfilename}','a') as sumfile:
            sumfile.write(f'IRC calculations submitted for {gfn}\n')
    # Launch conceptual DFT calculations
    elif ac == 'cdft':
        for fn in filenames[scope]:
            outfile = ro.Output(fn) 
            for i in range(1, 4): 
                jobname = f'{ac}_{outfile.base_name}_{i}'
                inpname = f'{jobname}.inp'
                chcdft = int(paramfile.charge) + i - 2
                multcdft = int(paramfile.mult) + i % 2
                #fm.writeg16inp(inpname,paramfile,12,outfile.coordinates,
                #               outfile.atnums,kw,chcdft,multcdft,ac=a)
                fm.writeinp(inpname,kw,outfile.coordinates,outfile.atnums,
                            paramfile,method,basis,chcdft,multcdft,
                            paramfile.cpus,solvent,disp,soft=software)
                fm.writecalcsubscript(paramfile,jobname,clust,soft=software)
                if not paramfile.woc:
                    os.system(f"{clust.subco} {jobname}.sub")
        with open(f'../{sumfilename}','a') as sumfile:
            sumfile.write(f'CDFT calculations submitted for {gfn}\n')
    # Launch NBO calculation
    elif ac == 'nbo':
        for fn in filenames[scope]:
            outfile = ro.Output(fn) 
            jobname = f'{ac}_{outfile.base_name}'
            inpname = f'{jobname}.inp'
            #fm.writeg16inp(inpname,paramfile,12,outfile.coordinates,
            #               outfile.atnums,kw,ac=True)
            fm.writeinp(inpname,kw,outfile.coordinates,outfile.atnums,
                        paramfile,method,basis,chcdft,multcdft,
                        paramfile.cpus,solvent,disp,soft=software)
            fm.writecalcsubscript(paramfile,jobname,clust,soft=software)
            if not paramfile.woc:
                os.system(f"{clust.subco} {jobname}g16.sub")
        with open(f'../{sumfilename}','a') as sumfile:
            sumfile.write(f'NBO calculations submitted for {gfn}\n')

# Call wheels.py
if paramfile.woc:
    dep = False
else:
    dep = True
scriptname = fm.writesubscript('wheels',paramfile,clust,dependency=dep)
os.system(f'{clust.subco} {scriptname}')


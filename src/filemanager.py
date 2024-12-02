import glob
import sys
import os
import re
import subprocess
import numpy as np
import spyrmsd
from spyrmsd import io, rmsd
fastcar_dir = "/home/hgeindre/ProjectFASTCAR/FASTCAR2"
import readoutput as ro

def writeconstraints(paramfile):
    """Write file constraints.inp with the constraints to be used in the
    CREST calculation"""
    with open('constraints.inp', 'w') as file:
        file.write('$constrain \n')
        file.write(f'  atoms: {paramfile.excluded_numbers_str} \n')
        file.write(f'  force constant= {paramfile.force} \n')
        file.write('  reference=struc.xyz \n')
        for bonds in paramfile.bonds_to_freeze:
            file.write(f'  distance: {bonds[0]}, {bonds[1]}, {bonds[2]} \n')
        for angles in paramfile.angles_to_freeze:
            file.write(f'  angle: {angles[0]}, {angles[1]}, {angles[2]}, '
                       f'{angles[3]} \n')
        for dihedral_angles in paramfile.dihedrals_to_freeze:
            file.write(f'  dihedral: {dihedral_angles[0]},'
                       f' {dihedral_angles[1]}, {dihedral_angles[2]}, '
                       f'{dihedral_angles[3]}, {dihedral_angles[4]} \n')
        file.write('$metadyn\n')
        file.write(f'  atoms: {paramfile.included_numbers_str} \n')
        file.write('$end') 

def writeparam(startoutput,paramfile):
    """Write a new parameter file parameters.tmp with the content of
    parameters.txt + experience name, value of first frequency charge and 
    multiplicity"""
    emp = ""
    heading = "General parameters"
    if startoutput.ts:
        ffreq = int(startoutput.freqs[0])
    else:
        ffreq = "none"
    if paramfile.check_st:
        exp_name = ''.join(''.join(startoutput.base_name.split('_')[1:])\
                   .split('-')[:-1])
        loop = re.search(r'Loop number (-?\d+)', paramfile.file_content)
        ln = loop.group(1)
        new_ln = str(int(ln) + 1)
        new_loop = loop.group().replace(ln,new_ln)
        modified_content = paramfile.file_content.replace(loop.group(),new_loop)
    else:
        exp_name = startoutput.base_name
        modified_content = f'{emp:-^50}\n{heading:-^50}\n{emp:-^50}\n'\
                         + f'Experience name {exp_name}\n'\
                         + f'First frequency {ffreq}\n'\
                         + f'Charge {startoutput.charge}\n'\
                         + f'Multiplicity {startoutput.multiplicity}\n'\
                         + f'Loop number {paramfile.check_st}\n'\
                         + paramfile.file_content
    if startoutput.software == "same":
        pattern = re.compile(r'[-]+DFT parameters[-]+\n[-]+\n', re.MULTILINE)
        match = pattern.search(modified_content)
        modified_content = paramfile.file_content[:match.end()] + '\n'\
                         + f'Software {startoutput.software}'\
                         + paramfile.file_content[match.end():]
    
    with open('parameters.tmp', 'w') as file:
        file.write(modified_content)
     
def writecrestsubscript(startoutput,pf):
    """Write bash script to submit CREST calculation"""
    multiplicityCREST = startoutput.multiplicity - 1
    with open('script_CREST.sub', 'w') as file:
        file.write('#!/bin/sh\n')
        file.write(f'#SBATCH --job-name=CREST_{startoutput.base_name} \n')
        file.write('#SBATCH --nodes=1 \n')
        file.write('#SBATCH --ntasks=6 \n')
        file.write(f'#SBATCH --output=CREST_{startoutput.base_name}.logfile \n')
        file.write('#SBATCH --time=48:00:00 \n')
        file.write('\n')
        file.write('cd ${SLURM_SUBMIT_DIR} \n')
        file.write('\n')
        file.write('#Loading modules \n')
        file.write('module load intel/2023.2.1\n')
        #file.write('module load intel/18.0.5.274 \n')
        #file.write('module load libraries/intelmpi/2018.4.274 \n')
        file.write('\n')
        file.write('#Launching calculation \n')
        file.write('\n')
        # Constrained CREST followed by pruning using CREGEN
        if not pf.woc:
            if startoutput.ts:
                file.write(f'{pf.version_crest} struc.xyz --T 6 --ewin '
                           f'{pf.ewin_crest} --uhf {multiplicityCREST} '
                           f'--chrg {startoutput.charge} {pf.solvent_crest}'
                           f' {pf.nci} -cinp constraints.inp --subrmsd > '
                           f'CrestAnalysis.txt\n')
                file.write(f'{pf.version_crest} coord -cregen '
                            'crest_conformers.xyz -ewin 30\n')
            # Unconstrained CREST
            else:
                file.write(f'{pf.version_crest} struc.xyz --T 6 --ewin '
                           f'{pf.ewin_crest} --uhf {multiplicityCREST} '
                           f'--chrg {startoutput.charge} {pf.nci} > '
                           f'CrestAnalysis.txt\n')
        # Calling engine.py for Pruning CREST conformers using SpyRMSD and 
        # submitting DFT calculations
        file.write(f'{fastcar_dir}/engine.py\n')
        file.write('\n')

def writeg16subscript(pf,jobname,clust,g16C=False):
    """Write script to submit gaussian calculation"""
    scriptname = f'{jobname}g16.sub'
    with open(scriptname,'w') \
            as gsub:
        gsub.write('#!/bin/sh\n')
        gsub.write(f'#SBATCH --job-name={jobname}\n')
        gsub.write('#SBATCH --nodes=1\n')
        if pf.node:
            if pf.node.replace("-", "", 1).isdigit():
                gsub.write(f'#SBATCH --exclude=node[{pf.node}]\n')
            if pf.node.isdigit():
                gsub.write(f'#SBATCH --exclude=node{pf.node}\n')
        gsub.write('#SBATCH --ntasks=12\n')
        gsub.write(f'#SBATCH --output={jobname}.logfile\n')
        gsub.write('#SBATCH --time=48:00:00\n')
        gsub.write('\n')
        gsub.write('#Defining Gaussian Parameters\n')
        if g16C:
            gsub.write(f'export g16root={clust.g16C_path}\n')
        else:
            gsub.write(f'export g16root={clust.g16_path}\n')
        gsub.write('source $g16root/bsd/g16.profile \n')
        gsub.write('export GAUSS_EXEDIR=$g16root\n')
        gsub.write('export LD_LIBRARY_PATH=$g16root/\n')
        gsub.write('export GAUSS_SCRDIR=${SLUMR_SUBMIT_DIR}\n')
        gsub.write('cd ${SLURM_SUBMIT_DIR}\n')
        gsub.write('\n')
        gsub.write('#Loading modules\n')
        #gsub.write('module load intel/18.0.5.274\n')
        #gsub.write('module load libraries/intelmpi/2018.4.274\n')
        gsub.write('module load intel/2023.2.1\n')
        gsub.write('export PATH=/opt/ohpc/pub/software/g16/:$PATH\n')
        gsub.write('export LD_LIBRARY_PATH=/opt/ohpc/pub/software/g16\n')
        gsub.write('\n')
        gsub.write('#Launching calculation\n')
        gsub.write(f'time g16 < {jobname}.inp > {jobname}.out\n')
        gsub.write('\n')

def writeorcasubscript(pf,scriptname,jobname,clust,g16C=False):
    """Write script to submit gaussian calculation"""
    with open(scriptname,'w') \
            as gsub:
        gsub.write('#!/bin/sh\n')
        gsub.write(f'#SBATCH --job-name={jobname}\n')
        gsub.write('#SBATCH --nodes=1\n')
        if pf.node:
            if pf.node.replace("-", "", 1).isdigit():
                gsub.write(f'#SBATCH --exclude=node[{pf.node}]\n')
            if pf.node.isdigit():
                gsub.write(f'#SBATCH --exclude=node{pf.node}\n')
        gsub.write('#SBATCH --ntasks=12\n')
        gsub.write(f'#SBATCH --output={jobname}.logfile\n')
        gsub.write('#SBATCH --time=48:00:00\n')
        gsub.write('\n')
        gsub.write('#Defining Gaussian Parameters\n')
        if g16C:
            gsub.write(f'export g16root={clust.g16C_path}\n')
        else:
            gsub.write(f'export g16root={clust.g16_path}\n')
        gsub.write('source $g16root/bsd/g16.profile \n')
        gsub.write('export GAUSS_EXEDIR=$g16root\n')
        gsub.write('export LD_LIBRARY_PATH=$g16root/\n')
        gsub.write('export GAUSS_SCRDIR=${SLUMR_SUBMIT_DIR}\n')
        gsub.write('cd ${SLURM_SUBMIT_DIR}\n')
        gsub.write('\n')
        gsub.write('#Loading modules\n')
        gsub.write('module load intel/2023.2.1\n')
        #gsub.write('module load intel/18.0.5.274\n')
        #gsub.write('module load libraries/intelmpi/2018.4.274\n')
        gsub.write('\n')
        gsub.write('#Launching calculation\n')
        gsub.write(f'time g16 < {jobname}.inp > {jobname}.out\n')
        gsub.write('\n')

def writecheckscript(pf):
    """Write script to call transmission script once all calculations have 
    ended"""
    name_script = 'script_check.sub'
    job_name = f'check_{pf.experience_number}'
    with open(f'{name_script}', 'w') as file:
        file.write('#!/bin/sh\n')
        file.write(f'#SBATCH --job-name={job_name}\n')
        file.write('#SBATCH --nodes=1 \n')
        file.write("#SBATCH --ntasks=1 \n")
        file.write(f'#SBATCH --output={job_name}.logfile \n')
        file.write('#SBATCH --time=48:00:00 \n')
        file.write('\n')
        file.write('cd ${SLURM_SUBMIT_DIR} \n')
        file.write('\n')
        file.write('#Loading modules \n')
        file.write('module load intel/2023.2.1\n')
        #file.write('module load intel/18.0.5.274 \n')
        #file.write('module load libraries/intelmpi/2018.4.274 \n')
        file.write('\n')
        file.write('#Launching calculation \n')
        file.write('\n')
        file.write(f'JobID=$(squeue --format="%.18i %.9P %.60j %.8u %.8T %.10M'
                   f' %.9l %.6D %R" -u $USER | grep "{pf.base_name_1}_'
                   f'{pf.experience_number}" | awk \'{{print $1}}\' | xargs | '
                   f'sed "s/ /,/g")\n')
        if pf.woc:
            file.write(f'sbatch {fastcar_dir}/transmission.py\n')
        else:
            file.write('sbatch --dependency=afterany:$JobID '
                      f'{fastcar_dir}/transmission.py\n')
        file.write('\n')

def writeabsscript(pf):
    """Write script to call abs"""
    name_script = 'script_abs.sub'
    job_name = f'abs_{pf.experience_number}'
    with open(f'{name_script}', 'w') as file:
       file.write('#!/bin/sh\n')
       file.write(f'#SBATCH --job-name={job_name}\n')
       file.write('#SBATCH --nodes=1 \n')
       file.write("#SBATCH --ntasks=1 \n")
       file.write(f'#SBATCH --output={job_name}.logfile \n')
       file.write('#SBATCH --time=48:00:00 \n')
       file.write('\n')
       file.write('cd ${SLURM_SUBMIT_DIR} \n')
       file.write('\n')
       file.write('#Loading modules \n')
       file.write('module load intel/2023.2.1\n')
       #file.write('module load intel/18.0.5.274 \n')
       #file.write('module load libraries/intelmpi/2018.4.274 \n')
       file.write('\n')
       file.write('#Launching calculation \n')
       file.write('\n')
       file.write(f'sbatch {fastcar_dir}/abs.py\n')
       file.write('\n')

def writeloopscript(pf):
    """Write script to call fastcar to loop for a check calculation"""
    name_script = 'script_loop.sub'
    nextloop = pf.check_st + 1
    job_name = f'{pf.experience_number}_loop{nextloop}'
    with open(f'{name_script}', 'w') as file:
       file.write('#!/bin/sh\n')
       file.write(f'#SBATCH --job-name={job_name}\n')
       file.write('#SBATCH --nodes=1 \n')
       file.write("#SBATCH --ntasks=1 \n")
       file.write(f'#SBATCH --output={job_name}.logfile \n')
       file.write('#SBATCH --time=48:00:00 \n')
       file.write('\n')
       file.write('cd ${SLURM_SUBMIT_DIR} \n')
       file.write('\n')
       file.write('#Loading modules \n')
       file.write('module load intel/2023.2.1\n')
       #file.write('module load intel/18.0.5.274 \n')
       #file.write('module load libraries/intelmpi/2018.4.274 \n')
       file.write('\n')
       file.write('#Launching calculation \n')
       file.write('\n')
       file.write(f'sbatch {fastcar_dir}/sparkplug.py\n')
       file.write('\n')

def writefinalscript(pf):
    """Call wheels.py once all calculations have ended"""
    with open('script_wheels.sub', 'w') as file:
        file.write('#!/bin/sh\n')
        file.write(f'#SBATCH --job-name=wheels_{pf.experience_number} \n')
        file.write('#SBATCH --nodes=1 \n')
        file.write("#SBATCH --ntasks=1 \n")
        file.write(f'#SBATCH --output=wheels_{pf.experience_number}'
                    '.logfile \n')
        file.write('#SBATCH --time=48:00:00 \n')
        file.write('\n')
        file.write('cd ${SLURM_SUBMIT_DIR} \n')
        file.write('\n')
        file.write('#Loading modules \n')
        file.write('module load intel/2023.2.1\n')
        #file.write('module load intel/18.0.5.274 \n')
        #file.write('module load libraries/intelmpi/2018.4.274 \n')
        file.write('\n')
        file.write('#Launching calculation \n')
        file.write('\n')
        if pf.ref_freq and pf.choice == 'irc':
            file.write('JobID_AddCalc=$(squeue --format="%.18i %.9P %.60j '
                       '%.8u %.8T %.10M %.9l %.6D %R" -u $USER | grep '
                      f'"{pf.choice}_{pf.experience_number}" | awk '
                       '\'{{print $1}}\' | xargs | sed "s/ /,/g")\n')
            if pf.woc:
                file.write(f'sbatch {fastcar_dir}/wheels.py')
            else:
                file.write('sbatch --dependency=afterany:$JobID_AddCalc '
                          f'{fastcar_dir}/wheels.py\n')
        #elif pf.choice == 'none':
        #    file.write(f'sbatch {fastcar_dir}/wheels.py')
        elif pf.choice != 'none':
            file.write('JobID_Freq=$(squeue --format="%.18i %.9P %.60j %.8u'
                       ' %.8T %.10M %.9l %.6D %R" -u $USER | grep '
                      f'\"freq_{pf.base_name_1}_{pf.experience_number}\" | awk '
                       '\'{{print $1}}\' | xargs | sed \"s/ /,/g\")\n')
            file.write('JobID_AddCalc=$(squeue --format="%.18i %.9P %.60j '
                       '%.8u %.8T %.10M %.9l %.6D %R" -u $USER | grep '
                      f'"{pf.choice}_{pf.experience_number}" | awk '
                       '\'{{print $1}}\' | xargs | sed "s/ /,/g")\n')
            if pf.woc:
                file.write(f'sbatch {fastcar_dir}/wheels.py')
            else:
                file.write('sbatch --dependency=afterany:$JobID_Freq,'
                          f'$JobID_AddCalc {fastcar_dir}/wheels.py\n')
        else:
            file.write('JobID_Freq=$(squeue --format="%.18i %.9P %.60j %.8u '
                       '%.8T %.10M %.9l %.6D %R" -u $USER | grep '
                      f'"freq_{pf.base_name_1}_{pf.experience_number}" | awk '
                       '\'{{print $1}}\' | xargs | sed "s/ /,/g")\n')
            if pf.woc:
                file.write(f'sbatch {fastcar_dir}/wheels.py')
            else:
                file.write('sbatch --dependency=afterany:$JobID_Freq '
                          f'{fastcar_dir}/wheels.py\n')
        file.write('\n')
    

def writeg16inp(inpname,pf,nprocs,geo,atnums,keywords,ch="def",mult="def",
                chkname="None",link1=False,keywords2="None",refine=False):
    """Write input for Gaussian"""
    if refine:
        base = pf.base_2
    else:
        base = pf.base_1
    fout = open(inpname, "w+")
    fout.write(f"%nprocshared={nprocs}\n")
    if chkname != "None":
        fout.write(f'%Chk={chkname}\n')
    fout.write("%mem=5GB\n")
    fout.write(f"#{keywords} {pf.functional} {base} {pf.dispersion}"
                    f" {pf.solvent}\n")
    fout.write("\n")
    fout.write(f"H2\n")
    fout.write(f"\n")
    if ch == "def":
        ch = pf.charge
    if mult == "def":
        mult = pf.multiplicity
    fout.write(f"{ch} {mult}\n")
    for i,item in enumerate(geo):
        fout.write(f'{atnums[i]} ')
        fout.write(' '.join(map(str, item)))
        fout.write('\n')
    fout.write('\n')
    if link1:
        fout.write('--Link1--\n')
        fout.write('%nprocshared=12\n')
        if chkname != "None":
            fout.write(f'%Chk={chkname}\n')
        fout.write(f'# {keywords2} {pf.functional} {base} {pf.dispersion}'
                   f' {pf.solvent}\n')
        fout.write('\n')
        fout.write('H2\n')
        fout.write('\n')
        fout.write(f'{pf.charge} {pf.multiplicity}\n')
        fout.write('\n')

def writeorcainp(inpname,pf,geo,keywords,ch="def",mult="def",refine=False):
    """Write input for ORCA"""
    if refine:
        base = pf.base_2
    else:
        base = pf.base_1
    fout = open(inpname, "w+")
    fout.write(f"! {keywords} {pf.functional} {pf.dispersion}"
                    f" {pf.solvent}\n")
    fout.write("! PrintBasis {base}\n")
    if 'optTS' in keywords:
        fout.write("%geom\n    Calc_Hess true\nend\n")
    fout.write("%output\n    print[p_mos] 1\nend\n")
    if ch == "def":
        ch = pf.charge
    if mult == "def":
        mult = pf.multiplicity
    fout.write(f"* xyz {ch} {mult}\n")
    for item in geo:
        fout.write(' '.join(map(str, item)))
        fout.write('\n')
    fout.write('*')
    #if link1:
    #    fout.write('--Link1--\n')
    #    fout.write('%nprocshared=12\n')
    #    if chk != "None":
    #        fout.write(f'%Chk={chkname}\n')
    #    fout.write(f'# {keywords2} {pf.functional} {base} {pf.dispersion}'
    #               f' {pf.solvent}\n')
    #    fout.write('\n')
    #    fout.write('H2\n')
    #    fout.write('\n')
        
def writegeolist(filename,geos,wm):
    with open(filename,wm) as gl:
        for fn in geos:
            gl.write(f'{len(geos[fn][1])}\n')
            gl.write(f'{fn}\n')
            for i,item in enumerate(geos[fn][1]):
                gl.write(f'{geos[fn][0][i]} ')
                gl.write(' '.join(map(str, item)))
                gl.write('\n')

def outstoxyzfile(filenames,xyzname):
    with open(xyzname,'w') as xyzfile:
        for fn in filenames:
            out = ro.Output(fn)
            xyzfile.write(str(out.num_atom))
            xyzfile.write('\n')
            xyzfile.write(fn)
            xyzfile.write('\n')
            for i,item in enumerate(out.coordinates):
                xyzfile.write(f'{out.atnums[i]} ')
                xyzfile.write(' '.join(map(str, item)))
                xyzfile.write('\n')

def writeresults(Energies,NPA,local_E,NBO):
    """Write results.txt with energies and potentially electrophilicity and
    bond energies"""
    with open('../results.txt', 'w') as log:
        log.write('--------------------------------------------------\n')
        log.write('----------Energies of final calculations----------\n')
        log.write('--------------------------------------------------\n')
        log.write('\n')
    
        for item in Energies:
            for element in item:
                if element:
                    log.write(str(element) + '\n')
            log.write('\n')
    
        if len(NPA) > 0:
            log.write('------------------------NPA------------------------')
            log.write('\n')
            log.write('\n')
            log.write('Atom   Number   Local electrophilicity\n\n')
            for item in local_E:
                for element in item:
                    log.write(str(element) + '        ')
                log.write('\n')
            log.write('\n')
    
        if len(NBO) > 0:
            log.write('------------------------NBO------------------------')
            log.write('\n')
            log.write('\n')
            log.write('Bond        Energie\n\n')
            for item in NBO:
                for element in item:
                    log.write(str(element) + '        ')
                log.write('\n')
            log.write('\n')

def writelogmin(E,garbageGeo,wm):
    """Add informations about DFT optimised minima in summary.log"""
    emp = ''
    if wm == 'w':
        with open('../structures.log', wm) as log:
            heading = 'Large base'
            log.write(f'{emp:-^50}\n')
            log.write(f'{heading:-^50}\n')
            log.write(f'{emp:-^50}\n\n\n')
            log.write('\n')
            heading = 'Min unique'
            log.write(f'{heading:-^50}\n')
            log.write('\n')
            for item in E:
                log.write(' '.join(map(str, item)))
                log.write('\n')
            log.write('\n')
            heading = 'Garbage'
            log.write(f'{heading:-^50}\n')
            for item in garbageGeo:
                log.write(' '.join(map(str, item)))
                log.write('\n')
    elif wm == 'a':
        pattern = re.compile(r'[-]+Garbage[-]+\n', re.MULTILINE)
        new_min_content = ''
        new_gar_content = ''
        for item in E:
            new_min_content += (' '.join(map(str, item)))
            new_min_content += ('\n')
        if len(garbageGeo) > 0:
            for item in garbageGeo:
                new_gar_content += (' '.join(map(str, item)))
                new_gar_content += ('\n')
        with open('../structures.log', 'r') as log:
            log_content = log.read()
            match = pattern.search(log_content)
            modified_content = log_content[:match.start()]\
                             + new_min_content + '\n'\
                             + log_content[match.start():]\
                             + new_gar_content
        with open('../structures.log', 'w') as log:
            log.write(modified_content)
            

def writelogTS(uniqueTS,doublonTS,garbageTS,failedTS,notcompletedTS,wm):
    """Add informations about DFT optimised transition states in summary.log"""
    emp = ''
    if wm == 'w':
        with open('../structures.log', wm) as log:
            log.write('\n')
            heading = 'Large base'
            log.write(f'{emp:-^50}\n')
            log.write(f'{heading:-^50}\n')
            log.write(f'{emp:-^50}\n\n\n')
            log.write('\n')
            heading = 'TS unique'
            log.write(f'{heading:-^50}\n')
            for item in uniqueTS:
                log.write(' '.join(map(str, item)))
                log.write('\n')
            log.write('\n')
            heading = 'Doublon'
            log.write(f'{heading:-^50}\n')
            if len(doublonTS) > 0:
                for item in doublonTS:
                    #log.write(' '.join(map(str, item)))
                    log.write(item)
                    log.write('\n')
            log.write('\n')
            heading = 'Garbage'
            log.write(f'{heading:-^50}\n')
            for item in garbageTS:
                log.write(' '.join(map(str, item)))
                log.write('\n')
            for item in failedTS:
                log.write(' '.join(map(str, item)))
                log.write('\n')
            for item in notcompletedTS:
                log.write(' '.join(map(str, item)))
                log.write('\n')
    elif wm == 'a':
        pattern_g = re.compile(r'[-]+Garbage[-]+\n', re.MULTILINE)
        pattern_d = re.compile(r'[-]+Doublon[-]+\n', re.MULTILINE)
        new_TS_content = ''
        new_db_content = ''
        new_gb_content = ''
        for item in uniqueTS:
            new_TS_content += (' '.join(map(str, item)))
            new_TS_content += ('\n')
        if len(doublonTS) > 0:
            for item in doublonTS:
                new_db_content += (' '.join(map(str, item)))
                new_db_content += ('\n')
        if len(garbageTS) > 0:
            for item in garbageTS:
                new_gb_content += (' '.join(map(str, item)))
                new_gb_content += ('\n')
        if len(failedTS) > 0:
            for item in failedTS:
                new_gb_content += (' '.join(map(str, item)))
                new_gb_content += ('\n')
        if len(notcompletedTS) > 0:
            for item in notcompletedTS:
                new_gb_content += (' '.join(map(str, item)))
                new_gb_content += ('\n')
        with open('../structures.log', 'r') as log:
            log_content = log.read()
            match_d = pattern_d.search(log_content)
            match_g = pattern_g.search(log_content)
            modified_content=log_content[:match_d.start()]\
                            +new_TS_content + '\n'\
                            +log_content[match_d.start():match_g.start()]\
                            +new_db_content + '\n'\
                            +log_content[match_g.start():]\
                            +new_gb_content 
        with open('../structures.log', 'w') as log:
            log.write(modified_content)
    
def crest_normterm(crestfile,ref_freq):
    """Check if crest calculation terminated normally"""
    with open(crestfile, 'r') as readfile:
        lines = readfile.readlines()
        last_line = lines[-1].rstrip()
        last_line_2 = lines[-2].rstrip()
    if ' CREST terminated normally.' in last_line:
        crestgood = True
        noreftopo = False
    elif 'or by using a method with fixed topology (GFN-FF).'  in last_line_2:
        crestgood = False
        noreftopo = True
    else:
        crestgood = False
        noreftopo = False
    return crestgood,noreftopo
        
def xyzspliter(xyzname,boutname="xyzfilenum"):
    """Function to split a concatenated XYZ file into separated XYZ files"""
    with open(xyzname, 'r') as rfile:
        lines = rfile.readlines()

    natoms = int(lines[0])
    ngeoms = len(lines) // (natoms + 2)

    for j in range(ngeoms):
        outname = f"{boutname}{j+1:04d}.xyz"
        with open(outname, "w") as ow:
            ow.write(str(natoms) + "\n \n")
            ow.write(lines[(j*(natoms + 2) + 2):((j + 1)*(natoms + 2))])

def xyzfinder(xyzname,n):
    """Function to read coordinate of a structure from an XYZ file containing
    several structures"""
    with open(xyzname, 'r') as rfile:
        lines = rfile.readlines()

    natoms = int(lines[0])
    return(lines[(n*(natoms + 2) + 2):((n + 1)*(natoms + 2))])

def xyzlistcleaner(filelist,RMSD_threshold):
    """Function identifying duplicate geometries in a list of XYZ files, 
       up to cutoff in RMSD compare to the first geometry only"""
    toremovefromlist = []
    RMSD_value = []
    rmsds = []
    if len(filelist) > 1:
        result = subprocess.check_output(f'python -m spyrmsd -m {filelist[0]} '
                                           'xyzfile*num*', shell=True)
        RMSD_value = result.split()
        RMSD_value = [item.decode() for item in RMSD_value]
        RMSD_value.extend(RMSD_value)

        for i in range(1,len(filelist)):
            if 0 < float(RMSD_value[i]) <= float(RMSD_threshold):
                toremovefromlist.append(filelist[i])
                rmsds.append(RMSD_value[i])

    return toremovefromlist,rmsds

def xyzlistcleanerbis(filename,RMSD_threshold,compstr):
    """Function identifying duplicate geometries in a list of XYZ files, 
       up to cutoff in RMSD compare to the first geometry only"""
    toremovefromlist = []
    RMSD_value = []
    rmsds = []
    compfiles = glob.glob(compstr)
    compfiles = sorted(compfiles)
    result = subprocess.check_output(f'python -m spyrmsd -m {filename} '
                                     f'{compstr}', shell=True)
    RMSD_value = result.split()
    RMSD_value = [item.decode() for item in RMSD_value]
    #RMSD_value.extend(RMSD_value)

    doublon = False
    mini = 999.0
    idfile = 'none'
    for f,cfile in enumerate(compfiles):
        if float(RMSD_value[f]) <= float(RMSD_threshold) and cfile != filename:
            doublon = True
            if float(RMSD_value[f]) < mini:
                mini = float(RMSD_value[f])
                idfile = cfile

    return doublon,idfile,mini

def geocomp(comp_mols,ref_mols,rmsd_thresh,strip=True):
    doublons = dict()
    nstrucs = len(comp_mols)
    if strip:
        for rmol in ref_mols:
            rmol.strip()
    rcoords = [np.array(mol.coordinates) for mol in ref_mols]
    ranum = ref_mols[0].atomicnums
    radj = ref_mols[0].adjacency_matrix
    for c,cmol in enumerate(comp_mols):
        if strip:
            cmol.strip()
        ccoords = cmol.coordinates
        canum = cmol.atomicnums
        cadj = cmol.adjacency_matrix
        rmsd_values = rmsd.symmrmsd(ccoords,rcoords,canum,ranum,cadj,radj, 
                                    minimize=True)
        mini = 999
        for r,rmsdv in enumerate(rmsd_values):
            if rmsdv < rmsd_thresh and rmsdv < mini:
                doublons[c] = (r,rmsdv)
                mini = rmsdv
    return doublons

def geocomp_onefile(comp_mols,rmsd_thresh,strip=True):
    doublons = dict()
    alreadyin = list()
    if strip:
        for cmol in comp_mols:
            cmol.strip()
    rcoords = [np.array(mol.coordinates) for mol in comp_mols]
    ranum = comp_mols[0].atomicnums
    radj = comp_mols[0].adjacency_matrix
    for c,cmol in enumerate(comp_mols):
        if c in alreadyin:
            continue
        ccoords = cmol.coordinates
        canum = cmol.atomicnums
        cadj = cmol.adjacency_matrix
        rmsd_values = rmsd.symmrmsd(ccoords,rcoords,canum,ranum,cadj,radj, 
                                    minimize=True)
        for r,rmsdv in enumerate(rmsd_values):
            if r == c:
                continue
            if rmsdv < rmsd_thresh:
                alreadyin.append(r)
                if not r in doublons:
                    doublons[r] = (c,rmsdv)
                else:
                    doublons[c] = (r,rmsdv)
    return doublons
    
def geocomp_onefile_complicated(comp_mols,rmsd_thresh,strip=True):
    doublons = list()
    similar = dict()
    alreadyin = list()
    if strip:
        for cmol in comp_mols:
            cmol.strip()
    rcoords = [np.array(mol.coordinates) for mol in comp_mols]
    ranum = comp_mols[0].atomicnums
    radj = comp_mols[0].adjacency_matrix
    for c,cmol in enumerate(comp_mols):
        #if c in alreadyin:
        #    continue
        ccoords = cmol.coordinates
        canum = cmol.atomicnums
        cadj = cmol.adjacency_matrix
        rmsd_values = rmsd.symmrmsd(ccoords,rcoords,canum,ranum,cadj,radj, 
                                    minimize=True)
        #print(c,rmsd_values)
        mini = 999
        idstruc = 'none'
        newlist1 = True
        for r,rmsdv in enumerate(rmsd_values):
            if r == c:
                continue
            if rmsdv < rmsd_thresh:
                newlist1 = False
                if rmsdv < mini:
                    idstruc = (r,rmsdv)
                    similar[c] = (r,rmsdv)
                    mini = rmsdv
                newlist2 = True
                for slist in doublons:
                    if r in slist and c not in slist:
                        slist.append(c)
                        newlist2 = False
                    elif c in slist and r not in slist:
                        newlist2 = False
                        slist.append(r)
                    elif c in slist and r in slist:
                        newlist2 = False
                if newlist2:
                    doublons.append([c,r])
        if newlist1:
            doublons.append([c])
                        
    return doublons,similar

def readsummary(lfn):
    nrjs = dict()
    nrjs[lfn] = dict()
    read = False
    with open(lfn,'r') as f:
        filelines = f.readlines()
        for line in filelines:
            if not read and 'TS unique' in line or 'Min unique' in line:
                if 'TS unique' in line:
                    structype = 'TS'
                elif 'Min unique' in line:
                    structype = 'min'
                read = True
            elif read and '------------' in line:
                read = False
            elif read and not line.isspace():
                ofn = line.split()[0]
                if structype == 'min':
                    nrj = float(line.split()[1])
                elif structype == 'TS':
                    nrj = float(line.split()[4])
                    freq = float(line.split()[8])
                nrjs[lfn][ofn] = nrj
    return nrjs
    
def writecdftlog(cdft_gds,cdft_funcs,ofiles,classif='charge'):
    cdft_fnames = [(r'f$^{+}$','f+'),(r'f$^{-}$','f-'),('\u0394f','delta_f'),
                   (r's$^{+}$','s+'),(r's$^{-}$','s-'),('\u0394s','delta_s'),
                   (u'\u0394\u03C1$_{elec}$','delta_rho_elec'),
                   (u'\u0394\u03C1$_{nuc}$','delta_rho_nuc')]
    fns = []
    for of in ofiles:
        fns.append((of.charge,of.file_name))
    fns.sort(key=lambda k: k[0]) 
    with open('cdft.log', 'w') as log:
        heading = f'File: {fns[1][1]}'
        emp = ""
        log.write(f'{heading:-^99}\n')
        log.write(f'Ionisation potential = {cdft_gds[0]} au.\n'
                  f'Electronic affinity = {cdft_gds[1]} au.\n' 
                  f'mu = {cdft_gds[2]} au.\n' 
                  f'mu+ = {cdft_gds[3]} au.\n' 
                  f'mu- = {cdft_gds[4]} au.\n' 
                  f'eta = {cdft_gds[5]} au.\n' 
                  f'omega = {cdft_gds[6]} au.\n\n')
        if classif == 'func':
            for f,func in enumerate(cdft_fnames):
                log.write(f'{func[1]}\n')
                log.write(f'{emp:<3}')
                for cht in cdft_funcs[f]:
                    log.write(f'{cht:>12}')
                log.write('\n')
                nats = len(cdft_funcs[f][cht])
                for i in range(0,nats):
                    log.write(f'{ofiles[1].atnums[i]:<3}')
                    for cht in cdft_funcs[f]:
                        log.write(f'{cdft_funcs[f][cht][i]:>12.6f}')
                    log.write('\n')
                log.write('\n')
        elif classif == 'charge':
            for cht in cdft_funcs[0]:
                log.write(f'{cht}\n')
                log.write(f'{emp:<3}')
                for f,func in enumerate(cdft_fnames):
                    funcstr = func[1].replace('delta_','d')
                    log.write(f'{funcstr:>12}')
                log.write('\n')
                nats = len(cdft_funcs[f][cht])
                for i in range(0,nats):
                    log.write(f'{ofiles[1].atnums[i]:<3}')
                    for f,func in enumerate(cdft_fnames):
                        log.write(f'{cdft_funcs[f][cht][i]:>12.6f}')
                    log.write('\n')
                log.write('\n')
                
def writecdftcsv(cdft_gds,cdft_funcs,ofiles,classif='charge'):
    cdft_fnames = [(r'f$^{+}$','f+'),(r'f$^{-}$','f-'),('\u0394f','delta_f'),
                   (r's$^{+}$','s+'),(r's$^{-}$','s-'),('\u0394s','delta_s'),
                   (u'\u0394\u03C1$_{elec}$','delta_rho_elec'),
                   (u'\u0394\u03C1$_{nuc}$','delta_rho_nuc')]
    fns = []
    for of in ofiles:
        fns.append((of.charge,of.file_name))
    fns.sort(key=lambda k: k[0]) 
    with open('cdft.csv', 'w') as log:
        log.write(f'{fns[1][1]};\n')
        log.write(f'Ionisation potential (au);{cdft_gds[0]}\n'
                  f'Electronic affinity (au);{cdft_gds[1]}\n' 
                  f'mu (au);{cdft_gds[2]}\n' 
                  f'mu+ (au);{cdft_gds[3]}\n' 
                  f'mu- (au);{cdft_gds[4]}\n' 
                  f'eta (au);{cdft_gds[5]}\n' 
                  f'omega (au);{cdft_gds[6]}\n\n')
        if classif == 'func':
            for f,func in enumerate(cdft_fnames):
                log.write(f'{func[1]};')
                for cht in cdft_funcs[f]:
                    log.write(f'{cht};')
                log.write('\n')
                nats = len(cdft_funcs[f][cht])
                for i in range(0,nats):
                    log.write(f'{ofiles[1].atnums[i]};')
                    for cht in cdft_funcs[f]:
                        log.write(f'{cdft_funcs[f][cht][i]};')
                    log.write('\n')
                log.write('\n')
        elif classif == 'charge':
            for cht in cdft_funcs[0]:
                log.write(f'{cht};')
                for f,func in enumerate(cdft_fnames):
                    log.write(f'{func[1]};')
                log.write('\n')
                nats = len(cdft_funcs[f][cht])
                for i in range(0,nats):
                    log.write(f'{ofiles[1].atnums[i]};')
                    for f,func in enumerate(cdft_fnames):
                        log.write(f'{cdft_funcs[f][cht][i]};')
                    log.write('\n')
                log.write('\n')
                
def readepc(keys='number'):
    epcPath = fastcar_dir + "/elementPlotCharac.data"
    try:
        elemCharacsFile = open(epcPath, 'r')
    except FileNotFoundError:
        print("File " + epcPath + "  not found")
        sys.exit(1)
    ecs = dict()
    ecsLines = elemCharacsFile.readlines()
    if keys == 'number':
        for n, l in enumerate(ecsLines):
            ecs[n + 1] = l.split()
    elif keys == 'letter':
        for n, l in enumerate(ecsLines):
            ecs[l.split()[0]] = [n+1] + l.split()[1:]
    for elt in ecs:
        ecs[elt][-1] = int(ecs[elt][-1])
    return ecs
        
    

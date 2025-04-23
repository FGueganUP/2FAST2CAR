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
        if paramfile.bonds:
            for bonds in paramfile.bonds_to_freeze:
                file.write(f'  distance: {bonds[0]}, {bonds[1]}, {bonds[2]} \n')
        if paramfile.angles:
            for angles in paramfile.angles_to_freeze:
                file.write(f'  angle: {angles[0]}, {angles[1]}, {angles[2]}, '
                           f'{angles[3]} \n')
        if paramfile.dihedrals:
            for dihedral_angles in paramfile.dihedrals_to_freeze:
                file.write(f'  dihedral: {dihedral_angles[0]},'
                           f' {dihedral_angles[1]}, {dihedral_angles[2]}, '
                           f'{dihedral_angles[3]}, {dihedral_angles[4]} \n')
        file.write('$metadyn\n')
        file.write(f'  atoms: {paramfile.included_numbers_str} \n')
        file.write('$end') 

def writeparam(startoutput,paramfile,clust):
    """Write a new parameter file parameters.tmp with the content of
    parameters.txt + experience name, value of imaginary frequency, charge and 
    multiplicity"""
    emp = ""
    heading = "general parameters"
    if startoutput.ts:
        ffreq = int(startoutput.freqs[0])
    else:
        ffreq = "none"
    if paramfile.loopn:
        exp_name = ''.join(''.join(startoutput.base_name.split('_')[1:])\
                   .split('-')[:-1])
        loop = re.search(r'loop number (-?\d+)', paramfile.file_content)
        ln = loop.group(1)
        new_ln = str(int(ln) + 1)
        new_loop = loop.group().replace(ln,new_ln)
        modified_content = paramfile.file_content.replace(loop.group(),new_loop)
    else:
        exp_name = startoutput.base_name
        modified_content = f'{emp:-^50}\n{heading:-^50}\n{emp:-^50}\n'\
                         + f'experience name {exp_name}\n'\
                         + f'imaginary frequency {ffreq}\n'\
                         + f'charge {startoutput.charge}\n'\
                         + f'multiplicity {startoutput.mult}\n'\
                         + f'loop number {paramfile.loopn}\n\n'\
                         + paramfile.file_content
        new_soft = f'software {paramfile.software}'
        soft = re.search(r'software (-?\d+)', paramfile.file_content)
        if soft:
            modified_content = modified_content.replace(soft.group(),
                                                        new_soft)
        else:
            method = re.search(r'method (.+)', modified_content)
            modified_content = modified_content[:method.start()]  +\
                               f'software {paramfile.software}\n'\
                               + modified_content[method.start():]
    
    with open('parameters.tmp', 'w') as file:
        file.write(modified_content)

def updateparam(startoutput,pf,tmpf,clust):
    """Update parameters.tmp with details on additional calculations for case
    when only additional calculations are performed"""
    keywords = ['additional calculation', 'scope', 'additional method', 
                'additional basis set', 'additional software', 
                'additional solvent', 'additional dispersion']
    new_content = tmpf.file_content
    for kw in keywords:
        match = re.search(fr'^{kw} (.+)', pf.file_content, flags=re.MULTILINE)
        tmp_match = re.search(fr'^{kw} (.+)', tmpf.file_content, 
                              flags=re.MULTILINE)
        if match:
            if tmp_match:
                new_content = new_content.replace(tmp_match.group(),
                                                  match.group())
            else:
                new_content += f'{match.group()}\n'
        elif tmp_match:
                new_content = new_content.replace(tmp_match.group(),'')
    with open('parameters.tmp', 'w') as file:
        file.write(new_content)
            
def writecrestsubscript(pf,startoutput,jobname,clust,script):
    soft_path = clust.path['crest']
    mult = startoutput.mult - 1
    with open(f'{fastcar_dir}/Models/crest.sub','r') as model:
        model_content = model.read()
    model_content = model_content.replace('{jobname}', jobname)
    model_content = model_content.replace('{cpus}', '6')
    model_content = model_content.replace('{time}', pf.crestime)
    model_content = model_content.replace('{soft_path}', 
                                         f'{soft_path}{pf.crest_version}')
    model_content = model_content.replace('{ewin}', str(pf.ewin_crest))
    model_content = model_content.replace('{mult}', str(mult))
    model_content = model_content.replace('{charge}', str(startoutput.charge))
    if pf.solvent_crest == 'none':
        model_content = model_content.replace('{solvent}', '')
    else:
        model_content = model_content.replace('{solvent}', pf.solvent_crest)
    model_content = model_content.replace('{nci}', pf.nci)
    constr = re.search(r'([^\\]\{)constraints(.*?)([^\\]\})', model_content, 
                      re.DOTALL)
    if startoutput.ts == '':
        model_content = model_content.replace(constr.group(),'')
    else:
        constr_match = constr.group()[2:-1]
        constr_match = constr_match.replace('constraints ','')
        model_content = model_content.replace(constr.group(),f'{constr_match}')
    if script != 'none':
        model_content = model_content.replace('{script}', 
                                             f'{fastcar_dir}/{script}')
    else:
        model_content = model_content.replace('{script}', '')
         
    scriptname = f'{jobname}.sub'
    with open(scriptname,'w') as scriptfile:
        scriptfile.write(model_content)

def writecalcsubscript(pf,jobname,clust,time='pf',soft='pf'):
    if time == 'pf':
        time = pf.optime
    if soft == 'pf':
        soft = pf.software
    soft_path = clust.path[soft]
    with open(f'{fastcar_dir}/Models/{soft}.sub','r') as model:
        model_content = model.read()
    model_content = model_content.replace('{jobname}', jobname)
    model_content = model_content.replace('{cpus}', str(pf.cpus))
    model_content = model_content.replace('{time}', time)
    model_content = model_content.replace('{soft_path}', soft_path)
    nodex = re.search(r'([^\\]\{)nodex(.*?)([^\\]\})', model_content, 
                      re.DOTALL)
    if pf.node == 'none':
        model_content = model_content.replace(nodex.group(),'')
    else:
        node_match = nodex.group()[2:-1]
        node_match = node_match.replace('nodex ','')
        node_match = node_match.replace('[nodenumbs]',f'node[{pf.node}]')
        model_content = model_content.replace(nodex.group(),f'\n{node_match}\n')
         
    scriptname = f'{jobname}.sub'
    with open(scriptname,'w') as scriptfile:
        scriptfile.write(model_content)

def writesubscript(script,pf,clust,dependency=False,time='48:00:00'):
    jobname = f'{script}_{pf.calc_name}'
    with open(f'{fastcar_dir}/Models/script.sub','r') as model:
        model_content = model.read()
    model_content = model_content.replace('{jobname}', jobname)
    model_content = model_content.replace('{time}', time)
    if dependency:
        model_content = model_content.replace('{dependency}', 
             'JobID=$(squeue --format="%.18i %.9P %.60j %.8u %.8T %.10M %.9l '
            f'%.6D %R" -u $USER | grep "{pf.calc_name}" | awk \'{{print $1}}\''
             ' | xargs | sed "s/ /,/g")\n')
        model_content = model_content.replace('{script}', f'{clust.subco} '
                   f'--dependency=afterany:$JobID {fastcar_dir}/{script}.py\n')
    else:
        model_content = model_content.replace('{dependency}', '')
        model_content = model_content.replace('{script}', 
            f'{clust.subco} {fastcar_dir}/{script}.py\n')
    scriptname = f'{jobname}.sub'
    with open(scriptname,'w') as scriptfile:
        scriptfile.write(model_content)
    return scriptname

def writeinp(inpname,keywords,geo,atnums,pf,meth='pf',basis='pf',ch='pf',
             mult='pf',cpus='pf',solvent='pf',disp='pf',soft='pf',chk='nope',
             intcoo='nope',supcalc=False,keywords2='none'):
    if meth == 'pf':
        meth = pf.method
    if basis == 'pf':
        basis = pf.basis_set
    if ch == 'pf':
        ch = pf.charge
    if mult == 'pf':
        mult = pf.mult
    if cpus == 'pf':
        cpus = pf.cpus
    if solvent == 'pf':
        solvent = pf.solvent
    if disp == 'pf':
       disp = pf.dispersion 
    if soft == 'pf':
        soft = pf.software
    if os.path.isfile(f'{fastcar_dir}/Models/{soft}_{keywords[1]}.inp'):
        model_path = f'{fastcar_dir}/Models/{soft}_{keywords[1]}.inp'
    else:
        model_path = f'{fastcar_dir}/Models/{soft}.inp'
    with open(model_path,'r') as model:
        model_content = model.read()
    model_content = model_content.replace('{cpus}', str(cpus))
    if chk == 'nope':
        model_content = model_content.replace('{chk}\n', '')
    else:
        model_content = model_content.replace('{chk}', f'%chk={chk}')
    model_content = model_content.replace('{keywords}', keywords[0])
    model_content = model_content.replace('{method}', meth)
    model_content = model_content.replace('{basis-set}', basis)
    model_content = model_content.replace('{charge}', str(ch))
    model_content = model_content.replace('{mult}', str(mult))
    if solvent == 'none':
        model_content = model_content.replace('{solvent}', '')
    else:
        model_content = model_content.replace('{solvent}', solvent)
    if disp == 'none':
        model_content = model_content.replace('{dispersion}', '')
    else:
        model_content = model_content.replace('{dispersion}', disp)
    ic_mod = re.search(r'([^\\]\{)int coord(.*?)([^\\]\})', 
                           model_content, re.DOTALL)
    if ic_mod:
        if intcoo == 'nope':
            model_content = model_content.replace(ic_mod.group(),'')
        else:
            model_content = model_content.replace(ic_mod.group(),'')
            ic_content = ''
            for ic_type in intcoo:
                for ic in ic_type:
                    ics = ic.split()
                    if len(ics) == 2:
                        ct = 'B'
                    if len(ics) == 3:
                        ct = 'A'
                    if len(ics) == 4:
                        ct = 'D'
                    ic_match = ic_mod.group()[2:-1]
                    ic_match = ic_match.replace('int coord ','')
                    ic_match = ic_match.replace('[T]',ct)
                    ic_match = ic_match.replace('[A1]',ics[0])
                    ic_match = ic_match.replace('[A2]',ics[1])
                    if ct == 'A' or ct == 'D':
                        ic_match = ic_match.replace('[A3]',ics[2])
                    else:
                        ic_match = ic_match.replace(' [A3]','')
                    if ct == 'D':
                        ic_match = ic_match.replace('[A4]',ics[3])
                    else:
                        ic_match = ic_match.replace(' [A4]','')
                    ic_match = ic_match.replace(r'\{',r'{')
                    ic_match = ic_match.replace(r'\}',r'}')
                    ic_content += ic_match + '\n'
            model_content = model_content[:ic_mod.start()+1] + ic_content +\
                            model_content[ic_mod.start()+1:]
        
    geo_text = ''
    for at,(x,y,z) in enumerate(geo):
        if at == len(geo)-1:
            geo_text += f'{atnums[at]} {x:.6f} {y:.6f} {z:.6f}'
        else:
            geo_text += f'{atnums[at]} {x:.6f} {y:.6f} {z:.6f}\n'
    model_content = model_content.replace('{geo}', geo_text)
    if 'wfx' in keywords[0]:
        wfxname = inpname.split('.')[0] + '.wfx'
        model_content = model_content.replace('{wfx}', str(wfxname))
    else:
        model_content = model_content.replace('{wfx}\n', '')
    sup_calc = re.search(r'([^\\]\{)sup calc(.*?)([^\\]\})', model_content, 
                         re.DOTALL)
    if sup_calc:
        if supcalc:
            sc_match = sup_calc.group()[2:-1]
            sc_match = sc_match.replace('sup calc ','')
            sc_match = sc_match.replace('[cpus]',str(cpus))
            if chk == 'nope':
                sc_match = sc_match.replace('[chk]\n', '')
            else:
                sc_match = sc_match.replace('[chk]', f'%chk={chk}')
            sc_match = sc_match.replace('[keywords]', keywords2)
            sc_match = sc_match.replace('[method]', meth)
            sc_match = sc_match.replace('[basis-set]', basis)
            sc_match = sc_match.replace('[charge]', str(ch))
            sc_match = sc_match.replace('[mult]', str(mult))
            if solvent == 'none':
                sc_match = sc_match.replace('[solvent]', '')
            else:
                sc_match = sc_match.replace('[solvent]', solvent)
            if disp == 'none':
                sc_match = sc_match.replace('[dispersion]', '')
            else:
                sc_match = sc_match.replace('[dispersion]', disp)
            model_content = model_content.replace(sup_calc.group(),sc_match)
        else: 
            model_content = model_content.replace(sup_calc.group(),'')
    with open(inpname,'w') as inpfile:
        inpfile.write(model_content)

def writegeolist(filename,geos,wm):
    """Write list a text file with all the geometries given, format is: 
    number of atoms 
    name of the structure 
    cartesian coordinates"""
    with open(filename,wm) as gl:
        for fn in geos:
            gl.write(f'{len(geos[fn][1])}\n')
            gl.write(f'{fn}\n')
            for i,item in enumerate(geos[fn][1]):
                gl.write(f'{geos[fn][0][i]} ')
                item = [f'{numb:.6f}' for numb in item]
                gl.write(' '.join(item))
                #gl.write(' '.join(map(str, item)))
                gl.write('\n')

def outstoxyzfile(filenames,xyzname):
    """Write a xyzfile from outputfiles, format is
    number of atoms
    name of the outputfile
    cartesian coordinates"""
    with open(xyzname,'w') as xyzfile:
        for fn in filenames:
            out = ro.Output(fn)
            xyzfile.write(str(out.num_atom))
            xyzfile.write('\n')
            xyzfile.write(fn)
            xyzfile.write('\n')
            for i,item in enumerate(out.coordinates):
                xyzfile.write(f'{out.atnums[i]} ')
                item = [f'{numb:.6f}' for numb in item]
                xyzfile.write(' '.join(item))
                #xyzfile.write(' '.join(map(str, item)))
                xyzfile.write('\n')

def writexyzfile(geo,atnums,filename,
                 commentline='file written with writexyzfile function'):
    """Write a xyzfile with one geometry, format is
    number of atoms
    comment
    cartesian coordinates"""
    with open(filename, 'w') as xyzfile:
        xyzfile.write(f'{len(geo)}\n')
        xyzfile.write(f'{commentline}\n')
        for at,c in enumerate(geo):
            xyzfile.write(f'{atnums[at]} {c[0]} {c[1]} {c[2]}\n')
        
    

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

def writelogmin(uniqueMin,duplicateMin,noMin,failedMin,wm):
    """Add informations about optimised minima in summary.log"""
    emp = ''
    if wm == 'w':
        with open('../structures.log', wm) as log:
            heading = 'Min unique'
            log.write(f'{heading:-^50}\n')
            log.write('\n')
            for item in uniqueMin:
                log.write(f'{item[0]} E(SCF) = {item[1]}')
                #log.write(' '.join(map(str, item)))
                log.write('\n')
            log.write('\n')
            heading = 'Duplicates'
            log.write(f'{heading:-^50}\n')
            if len(duplicateMin) > 0:
                for item in duplicateMin:
                    log.write(f'{item[0]} E SCF = {item[1]} id to {item[2]} '
                              f'rmsd = {item[3]}')
                    log.write('\n')
            log.write('\n')
            heading = 'Garbage'
            log.write(f'{heading:-^50}\n')
            for item in noMin:
                log.write(f'{item[0]} E(SCF) = {item[1]}, First frequency = '
                          f'{item[2]}')
                log.write(' '.join(item))
                #log.write(' '.join(map(str, item)))
                log.write('\n')
            for item in failedMin:
                log.write(' '.join(item))
                #log.write(' '.join(map(str, item)))
                log.write('\n')
    elif wm == 'a':
        pattern_g = re.compile(r'[-]+Garbage[-]+\n', re.MULTILINE)
        pattern_d = re.compile(r'[-]+Duplicates[-]+\n', re.MULTILINE)
        new_min_content = ''
        new_db_content = ''
        new_gb_content = ''
        for item in uniqueMin:
            new_min_content += f'{item[0]} E(SCF) = {item[1]}'
            #new_min_content += (' '.join(map(str, item)))
            new_min_content += ('\n')
        if len(duplicateMin) > 0:
            for item in duplicateMin:
                new_db_content += (f'{item[0]} E SCF = {item[1]} id to '
                                   f'{item[2]} rmsd = {item[3]}')
                #new_db_content += (' '.join(map(str, item)))
                new_db_content += ('\n')
        if len(noMin) > 0:
            for item in noMin:
                new_gb_content += (f'{item[0]} E(SCF) = {item[1]}, First '
                                   f'frequency = {item[2]}')
                #new_gb_content += (' '.join(map(str, item)))
                new_gb_content += ('\n')
        if len(failedMin) > 0:
            for item in failedMin:
                new_gb_content += (' '.join(item))
                #new_gb_content += (' '.join(map(str, item)))
                new_gb_content += ('\n')
        with open('../structures.log', 'r') as log:
            log_content = log.read()
            match_d = pattern_d.search(log_content)
            match_g = pattern_g.search(log_content)
            modified_content=log_content[:match_d.start()]\
                            +new_min_content + '\n'\
                            +log_content[match_d.start():match_g.start()]\
                            +new_db_content + '\n'\
                            +log_content[match_g.start():]\
                            +new_gb_content 
        with open('../structures.log', 'w') as log:
            log.write(modified_content)
            

def writelogTS(uniqueTS,duplicateTS,garbageTS,noTS,failedTS,wm):
    """Add informations about optimised transition states in summary.log"""
    emp = ''
    if wm == 'w':
        with open('../structures.log', wm) as log:
            heading = 'TS unique'
            log.write(f'{heading:-^50}\n')
            for item in uniqueTS:
                log.write(f'{item[0]} E(SCF) = {item[1]}, imaginary frequency '
                          f'= {item[2]}\n')
            log.write('\n')
            heading = 'Duplicates'
            log.write(f'{heading:-^50}\n')
            if len(duplicateTS) > 0:
                for item in duplicateTS:
                    log.write(f'{item[0]} E(SCF) = {item[1]}, imaginary '
                              f'frequency = {item[2]} id to {item[3]}, rmsd = '
                              f'{item[4]}\n')
            log.write('\n')
            heading = 'Garbage'
            log.write(f'{heading:-^50}\n')
            for item in garbageTS:
                log.write(f'{item[0]} E(SCF) = {item[1]}, imaginary frequency '
                          f'= {item[2]}\n')
            for item in noTS:
                log.write(f'{item[0]} E(SCF) = {item[1]}, imaginary frequency '
                          f'= {item[2]}\n')
            for item in failedTS:
                log.write(f'{item[0]} {item[1]}\n')
    elif wm == 'a':
        pattern_g = re.compile(r'[-]+Garbage[-]+\n', re.MULTILINE)
        pattern_d = re.compile(r'[-]+Duplicates[-]+\n', re.MULTILINE)
        new_TS_content = ''
        new_db_content = ''
        new_gb_content = ''
        for item in uniqueTS:
            new_TS_content += (f'{item[0]} E(SCF) = {item[1]}, imaginary '
                               f'frequency = {item[2]}')
            #new_TS_content += (' '.join(map(str, item)))
            new_TS_content += ('\n')
        if len(duplicateTS) > 0:
            for item in duplicateTS:
                new_db_content += (f'{item[0]} E(SCF) = {item[1]}, first '
                                   f'frequency = {item[2]} id to {item[3]}, '
                                   f'rmsd = {item[4]}\n')
        if len(garbageTS) > 0:
            for item in garbageTS:
                new_gb_content += (f'{item[0]} E(SCF) = {item[1]}, '
                                   f'imaginary frequency = {item[2]}\n')
        if len(noTS) > 0:
            for item in noTS:
                new_gb_content += (f'{item[0]} E(SCF) = {item[1]}, '
                                   f'imaginary frequency = {item[2]}\n')
        if len(failedTS) > 0:
            for item in failedTS:
                new_gb_content += (' '.join(item))
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
    ecs = readepc(keys='letter')
    n -= 1
    with open(xyzname, 'r') as rfile:
        lines = rfile.readlines()

    natoms = int(lines[0])
    geo = list()
    atnums = list()
    for line in lines[(n*(natoms + 2) + 2):((n + 1)*(natoms + 2))]:
        line.strip('\n')
        geo.append([float(c) for c in line.split()[1:]])
        try:
            atnum = int(line.split()[0])
        except:
            atnum = ecs[line.split()[0]][0]
        atnums.append(atnum)
    return(geo,atnums)

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

    duplicate = False
    mini = 999.0
    idfile = 'none'
    for f,cfile in enumerate(compfiles):
        if float(RMSD_value[f]) <= float(RMSD_threshold) and cfile != filename:
            duplicate = True
            if float(RMSD_value[f]) < mini:
                mini = float(RMSD_value[f])
                idfile = cfile

    return duplicate,idfile,mini

def geocomp(comp_mols,ref_mols,rmsd_thresh,strip=True):
    """Compare two sets of geometries with spyrmsd, return a list of files
    that are identical to another"""
    duplicates = dict()
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
                duplicates[c] = (r,rmsdv)
                mini = rmsdv
    return duplicates

def geocomp_onefile(comp_mols,rmsd_thresh,strip=True):
    """Compare all geometries in a set, return list of structures that are 
    identical to another one"""
    duplicates = dict()
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
                if not r in duplicates:
                    duplicates[r] = (c,rmsdv)
                else:
                    duplicates[c] = (r,rmsdv)
    return duplicates
    
def geocomp_onefile_complicated(comp_mols,rmsd_thresh,strip=True):
    """Useless"""
    duplicates = list()
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
                for slist in duplicates:
                    if r in slist and c not in slist:
                        slist.append(c)
                        newlist2 = False
                    elif c in slist and r not in slist:
                        newlist2 = False
                        slist.append(r)
                    elif c in slist and r in slist:
                        newlist2 = False
                if newlist2:
                    duplicates.append([c,r])
        if newlist1:
            duplicates.append([c])
                        
    return duplicates,similar

def readsummary(lfn):
    """Read informations from structures.log (need to change that function 
    name)"""
    nrjs = dict()
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
                nrj = float(line.split()[3].strip(','))
                if structype == 'TS':
                    freq = float(line.split()[7])
                nrjs[ofn] = nrj
    return nrjs
    
def writecdftlog(cdft_gds,cdft_funcs,ofiles,classif='charge'):
    """Write text file with results of CDFT calculations"""
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
                  #f'mu+ = {cdft_gds[3]} au.\n' 
                  #f'mu- = {cdft_gds[4]} au.\n' 
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
    """Write csv file with results of CDFT calculations"""
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
    """Read elementplotCharac.data file containing informations on the elements
    of the periodic table, return a dict() with the informations either with
    the numbers of the elements as keys or with the symbols of the elements as 
    keys"""
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
        
def comptheolvl(output,paramfile):
    """Check if given outputfile as the same level of theory than the one that
    will be used for opt calculations, if yes the outputfile is included
    in the results, if not it is reoptimised"""
    reopt = False
    if not paramfile.method in output.method or \
       paramfile.basis_set != output.basis_set or \
       paramfile.solvent != output.solvent or \
       paramfile.software != output.software:
        print('Level of theory for input file seems different that the one '
              'specified in parameters.txt (be careful, dispersion is not '
              'monitored)')
        print(f'Input file: {output.method} / {output.basis_set}, '
              f'solvent = {output.solvent}, {output.software}')
        print(f'Parameters: {paramfile.method} / {paramfile.basis_set} '
              f'solvent = {paramfile.solvent}, {paramfile.software}')
        ans = input('Do you confirm it is different? (yes/no)\n')
        while True:
            if ans == 'yes' or ans == 'y':
                print('Thanks, input geometry will be reoptimised at '
                     f'{paramfile.method} / {paramfile.basis_set} level with '
                     f'{paramfile.software} and incorporated in results')
                reopt = True
                break
            elif ans == 'no' or ans == 'n':
                print('Thanks, input calculation will be incorporated in '
                      'results')
                reopt = False
                break
            else:
                ans = input('Please answer "yes" if you want the structure '
                            ' given at the start to be reoptimised at the '
                            'level of theory specified in the parameters file '
                            'or "no" if you want the output file given at the '
                            'start to be directly incorporated in the results\n'
                           )
    return reopt
    

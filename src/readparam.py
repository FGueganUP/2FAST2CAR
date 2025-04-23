import sys
import re
import misc
import math
import glob
fastcar_dir = "/home/hgeindre/ProjectFASTCAR/FASTCAR2"

class Parameters:
    """Class that is used to read and store the informations in the parameters
    text file"""

    def __init__(self, file_name,startoutput=None,clust=None):
        """Open the file and read """
        self.fn = file_name
        self.loopn = "unknown"
        if self.fn.split('.')[-1] == 'txt':
            with open(self.fn,'r') as file_brut:
                self.file_content = file_brut.read().lower()
            with open(self.fn,'r') as file_brut:
                self.file_lines = file_brut.readlines()
                self.file_lines = [l.lower() for l in self.file_lines]
            self.checkfile(startoutput,clust)
        elif self.fn.split('.')[-1] == 'tmp':
            with open(self.fn,'r') as file_brut:
                self.file_content = file_brut.read()
            with open(self.fn,'r') as file_brut:
                self.file_lines = file_brut.readlines()
                self.file_lines = [l for l in self.file_lines]
            self.getgeneral()
        else:
            print(f"File {self.fn} should not be named like that")
            sys.exit(1)

    def checkfile(self,startoutput,clust):
        """Check if all informations needed are in the parameters file and if
        given informations are in the correct format"""
        #compsofts = ['g16','orca']
        adcalcs = ['single-point','wfx','reopt','irc','nbo','cdft','none']
        mkwds = ['method', 'basis set']
        gkwds = ['max loops', 'ewin', 'bond constrained', 'angle constrained',
                 'dihedral angle constrained', 'force constant', 'nci', 
                 'crest version', 'crest solvent', 'rmsd threshold', 
                 'crest time limit', 'method', 'basis set', 'software', 
                 'dispersion', 'opt solvent', 'time limit', 'cpus', 
                 'nodes excluded', 'constrained opt', 'ts screening', 
                 'active atoms', 'additional calculation', 'scope', 
                 'additional method', 'additional basis set', 
                 'additional software', 'additional solvent', 
                 'additional dispersion']
        tsc_opt = ['activats', 'scalprod']
        
        # Search for mandatory keywords
        for kw in mkwds:
            res = re.search(f'{kw}(.+)', self.file_content)
            if res:
                pass
            else:
                print(f'Unable to find {kw} keyword in {self.fn}')
                sys.exit(1)   
            
        # Check if wrong keywords
        kwds = re.findall(r'^(?:\b[a-zA-Z]+\b *){1,3}(?: +|$)', 
                          self.file_content, re.MULTILINE)
        for kwd in kwds:
            kwd = kwd.strip()
            good = False
            for gkwd in gkwds:
                if kwd.startswith(gkwd):
                    good = True
                    break
            if not good:
                print(f'Keyword {kwd} is not correct')
                sys.exit(1)
                
        # Check if given keywords are correct
        loop = re.search(r'^max loops (.+)', self.file_content,
                         flags=re.MULTILINE)
        if loop:
            if loop.group(1).strip() == 'none':
                pass
            else:
                try:
                    loop = int(loop.group(1).strip())
                except ValueError:
                    print('Max loops should be none or an integer')
                    sys.exit(1)
                if loop < 0:
                    print('Max loops should not be negative')
                    sys.exit(1)

        crest_version = re.search(r'^crest version(.+)', self.file_content,
                                  flags=re.MULTILINE)
        crestvs = ['default','3.0','continous release']
        if crest_version:
            crest_version = crest_version.group(1).strip() 
            if crest_version not in crestvs:
                print('Please provide a suitable CREST version:')
                for v in crestvs:
                    print(f'-{v}')
                sys.exit(1)
        ewin_crest = re.search(r'^ewin(.+)', self.file_content,
                               flags=re.MULTILINE)
        if ewin_crest:
            try:
                ewin_crest = float(ewin_crest.group(1).strip())
            except ValueError:
                print('Ewin should be a float')
                sys.exit(1)
        solvent_crest = re.search(r'^crest solvent(.+)',self.file_content,
                                  flags=re.MULTILINE)
        if solvent_crest:
            solvent_crest = solvent_crest.group(1).strip()
            if re.search(r'--alpb', solvent_crest) or \
               re.search(r'--gbsa', solvent_crest) or \
               solvent_crest == 'none':
                pass
            else:
                print('Please provide a suitable CREST solvent.')
                sys.exit(1)
        nci = re.search(r'^nci(.*)', self.file_content,
                        flags=re.MULTILINE)
        if nci:
            if len(nci.group().split()) > 2:
                print('Too many arguments in NCI option.')
                sys.exit(1)
            elif len(nci.group().split()) == 2 and nci.group(1).strip() \
                                                   != 'none':
                try:
                    nci_value = float(nci.group(1).strip())
                except ValueError:
                    print('Please provide a suitable NCI scale factor.')
                    sys.exit(1)  
        bonds = re.findall('bond constrained(.+)', self.file_content)
        angles = re.findall(r'^angle constrained(.+)', self.file_content, 
                            re.MULTILINE)
        dihedrals = re.findall('dihedral angle constrained (.+)', 
                               self.file_content)
        if bonds:
            for bond in bonds:
                atoms = bond.split()
                if len(atoms) != 2:
                    print('Issue with bond constraint format')
                    sys.exit(1)
                try:
                    atom1 = int(atoms[0])
                    atom2 = int(atoms[1])
                except ValueError:
                    print('Atoms number for bond constraint should be ints.')
                    sys.exit(1)
        if angles:
            for angle in angles:
                atoms = angle.split()
                if len(atoms) != 3:
                    print('Issue with angle constraint format')
                    sys.exit(1)
                try:
                    atom1 = int(atoms[0])
                    atom2 = int(atoms[1])
                    atom3 = int(atoms[2])
                except ValueError:
                    print('Atoms number for angle constraint should be ints.')
                    sys.exit(1)
        if dihedrals:
            for dihedral in dihedrals:
                atoms = dihedral.split()
                if len(atoms) != 4:
                    print('Issue with dihedral angle constraint format')
                    sys.exit(1)
                try:
                    atom1 = int(atoms[0])
                    atom2 = int(atoms[1])
                    atom3 = int(atoms[2])
                    atom4 = int(atoms[3])
                except ValueError:
                    print('Atoms number for dihedral angle constraint should '
                          'be ints.')
                    sys.exit(1)
            
            force = re.search(r'^force constant(.+)', self.file_content,
                              flags=re.MULTILINE)
            if force:
                force = force.group(1).strip()
                try:
                    force_value = float(force)
                except ValueError:
                    print('Force constant should be a float')
                    sys.exit(1)
            else:
                print('Unable to find force constant information in the '
                      'parameters file.')
                sys.exit(1)      

        RMSD = re.search(r'^rmsd threshold(.+)' , self.file_content,
                         flags=re.MULTILINE)
        if RMSD:
            if len(RMSD.group(1).split()) == 1:
                try:
                    RMSD = float(RMSD.group(1))
                except ValueError:
                    print('RMSD threshold should be a float')
                    sys.exit(1)
            elif len(RMSD.group(1).split()) == 2: 
                try:
                    RMSD_crest = float(RMSD.group(1).split()[0])
                    RMSD_opt = float(RMSD.group(1).split()[1])
                except ValueError:
                    print('RMSD threshold should be a float')
                    sys.exit(1)
            else:
                print('Please provide one or two RMSD values.')
                sys.exit(1) 
        crest_time = re.search(r'^crest time limit (.+)', self.file_content,
                         flags=re.MULTILINE)
        if crest_time:
            crest_time = crest_time.group(1).strip()
            goodformat = misc.timeformat(crest_time)
            if not goodformat:
                print('Crest time limit format incorrect')
                sys.exit(1)

        cpus = re.search(r'^cpus (.+)', self.file_content, flags=re.MULTILINE)
        if cpus:
            try:
                cpus = int(cpus.group(1).strip())
            except ValueError:
                print('cpus should be an integer')
                sys.exit(1)

        time = re.search(r'^time limit (.+)', self.file_content,
                         flags=re.MULTILINE)
        if time:
            optime = time.group(1).strip()
            goodformat = misc.timeformat(optime)
            if not goodformat:
                print('Time limit format incorrect')
                sys.exit(1)

        tsc = re.search(r'^TS screening (.+)' , self.file_content,
                        flags=re.MULTILINE)
        if tsc:
            tsc = tsc.group(1).strip()
            if tsc not in tsc_opt:
                print(f'TS screening should be one of {tsc_opt}.')
                sys.exit(1)   

        actats = re.search(r'^active atoms(.+)' , self.file_content,
                           flags=re.MULTILINE)
        if actats:
            for elt in actats.group(1).split():
                try:
                    nat = int(elt)
                except ValueError:
                    print('Active atoms should be integers')
                    sys.exit(1)
        elif tsc:
            if tsc == 'activats':
                print('Active atoms have been chosen for TS screening so they '
                      'should be specified in the parameters file with the '
                      'keyword active atoms')
                sys.exit(1)
                

        adcalc = re.search(r'^additional calculation (.+)', self.file_content, 
                           flags=re.MULTILINE)
        if adcalc:
            adcalc = [ac.strip() for ac in adcalc.group(1).split()]
            for ac in adcalc: 
                if ac not in adcalcs:
                    print("Additional calculation should be one of the "
                          "following:\n")
                    for gac in adcalcs:
                        print(f"- {gac}")
                    sys.exit(1) 

        scope = re.search(r'^scope (.+)', self.file_content,flags=re.MULTILINE)
        if scope:
            #scope = scope.group(1).strip()
            scope = [sc.strip() for sc in scope.group(1).split()]
            for sc in scope:
                if sc != "all" and sc != "lowest":
                    print("Scope should be all or lowest\n")
                    sys.exit(1) 

        self.software = re.search(r'^software(.+)', self.file_content, 
                                  flags=re.MULTILINE)
        inp_models = [n.split('/')[-1].split('.')[0] for n in 
                      glob.glob(f'{fastcar_dir}/Models/*.inp')]
        sub_models = [n.split('/')[-1].split('.')[0] for n in 
                      glob.glob(f'{fastcar_dir}/Models/*.sub')]
        aliases = list()
        for alias in inp_models:
            if alias in sub_models:
                aliases.append(alias)
        if self.software:
            if self.software.group(1).strip() == 'same':
                for soft_options in clust.default:
                    if soft_options in startoutput.software.lower():
                        soft = clust.default[soft_options]
                        break
                else:
                    print(f"There is no default version for the software "
                           "{startoutput.software} in the file config.txt")
                    sys.exit(1)
                if soft not  in aliases:
                    print(f'Software is unknown, it must be one of {aliases}.')
                    sys.exit(1)   
            elif self.software.group(1).strip() not in aliases:
                print(f'Software is unknown, it must be one of {aliases}.')
                sys.exit(1)   
        else:
            self.software = 'same'
        self.ad_software = re.search(r'^additional software(.+)', 
                                    self.file_content, flags=re.MULTILINE)
        if self.ad_software:
            if self.ad_software.group(1).strip() == 'same':
                for soft_options in clust.default:
                    if soft_options in startoutput.software.lower():
                        soft = clust.default[soft_options]
                        break
                else:
                    print(f"There is no default version for the software "
                           "{startoutput.software} in the file config.txt")
                    sys.exit(1)
                if soft not  in aliases:
                    print(f'Software is unknown, it must be one of {aliases}.')
                    sys.exit(1)   
            elif self.ad_software.group(1).strip() not in aliases:
                print(f'Software is unknown, it must be one of {aliases}.')
                sys.exit(1)   

    def getgeneral(self):
        """Read general informations about calculation"""
        calc_name = re.search(r'^experience name (.+)',self.file_content,
                              flags=re.MULTILINE)
        self.calc_name = calc_name.group(1).strip()
    
        charge = re.search(r'^charge (-?\d+)', self.file_content,
                           flags=re.MULTILINE)
        self.charge = charge.group(1).strip()

        mult = re.search(r'^multiplicity (-?\d+)', self.file_content, 
                         flags=re.MULTILINE)
        self.mult = mult.group(1).strip()

        ref_freq = re.search(r'^imaginary frequency (.+)', self.file_content,
                             flags=re.MULTILINE)
        if ref_freq.group(1).strip() == 'none':
            self.ref_freq = ''
        else:
            self.ref_freq = int(ref_freq.group(1))

        loopn = re.search(r'^loop number (.+)', self.file_content, 
                          flags=re.MULTILINE)
        self.loopn = int(loopn.group(1).strip())


    def getloop(self):
        """Read information about loops"""
        # Write WoC True in parameters file to simulate a run of fastcar but 
        # without actually doing all the calculations. Calculations must have
        # already be done. It's a debug feature
        loop = re.search(r'^max loops (.+)', self.file_content,
                         flags=re.MULTILINE)
        if not loop:
            self.loop = 10
        else:
            self.loop = int(loop.group(1).strip())
        
    def getcrest(self):
        """Read: 
        - Which version of CREST to use, 
        - The size of the energy window, if none, take the default value
        6 kcal/mol
        - NCI
        """
        crest_version = re.search(r'^crest version(.+)', self.file_content,
                                  flags=re.MULTILINE)
        if crest_version:
            self.crest_version = crest_version.group(1).strip() 
            if self.crest_version == 'default': 
                self.crest_version = 'crest'
            elif self.crest_version == '3.0':
                self.crest_version = 'crest3'
            elif self.crest_version == 'continous release':
                self.crest_version = 'crest_continous_release'
        else:
            self.crest_version = 'crest'

        self.ewin_crest = re.search(r'^ewin(.+)', self.file_content,
                                    flags=re.MULTILINE)
        if self.ewin_crest:
            self.ewin_crest = float(self.ewin_crest.group(1).strip())
        else:
            self.ewin_crest = 6

        self.solvent_crest = re.search(r'^crest solvent(.+)',self.file_content,
                                       flags=re.MULTILINE)
        if self.solvent_crest:
            self.solvent_crest = self.solvent_crest.group(1).strip()
        else:
            self.solvent_crest = 'none'

        nci = re.search(r'^nci(.*)', self.file_content,
                        flags=re.MULTILINE)
        if nci:
            if len(nci.group().split()) == 1:
                self.nci = '--nci'
            elif len(nci.group().split()) == 2:
                if nci.group(1) == 'none':
                    self.nci = ''
                else: 
                    nci_value = float(nci.group(1))
                    self.nci = f'--nci --wscal {self.nci_value}'
        else:
            self.nci = ''

        crest_time = re.search(r'^crest time limit (.+)', self.file_content,
                         flags=re.MULTILINE)
        if not crest_time:
            self.crestime = '48:00:00'
        else:
            self.crestime = crest_time.group(1).strip()

    def getconstraints(self,startoutput):
        """Read constraints from parameters file"""
        self.bonds = re.findall('bond constrained(.+)', self.file_content)
        self.angles = re.findall(r'^ *angle constrained(.+)', self.file_content, 
                                 re.MULTILINE)
        self.dihedrals = re.findall('dihedral angle constrained (.+)', 
                                    self.file_content)
        
        excluded_numbers = []
        coords = startoutput.coordinates
        if self.bonds:
            self.bonds_to_freeze = []
            excluded_numbers_bonds =  []
            for bond in self.bonds:
                atoms = bond.split()
                atom1 = list(map(float, coords[int(atoms[0]) - 1]))
                atom2 = list(map(float, coords[int(atoms[1]) - 1]))
                distance = math.dist(atom1, atom2)
                self.bonds_to_freeze.append((int(atoms[0]), int(atoms[1]), 
                                             f"{distance:.3f}"))
            excluded_numbers += [item for sublist in self.bonds_to_freeze for \
                                 item in sublist[:2]]

        if self.angles:
            self.angles_to_freeze = []
            excluded_numbers_angle =  []
            for angle in self.angles:
                atoms = angle.split()
                atom1 = list(map(float, coords[int(atoms[0]) - 1]))
                atom2 = list(map(float, coords[int(atoms[1]) - 1]))
                atom3 = list(map(float, coords[int(atoms[2]) - 1]))
                angle = misc.calculate_angle(atom1, atom2, atom3)
                self.angles_to_freeze.append((int(atoms[0]),int(atoms[1]), 
                                              int(atoms[2]),f"{angle:.2f}"))
            excluded_numbers += [item for sublist in self.angles_to_freeze for \
                                 item in sublist[:3]]

        if self.dihedrals:
            self.dihedrals_to_freeze = []
            excluded_numbers_dihedral = []
            for dihedral in self.dihedrals:
                atoms = dihedral.split()
                atom1 = list(map(float, coords[int(atoms[0]) - 1]))
                atom2 = list(map(float, coords[int(atoms[1]) - 1]))
                atom3 = list(map(float, coords[int(atoms[2]) - 1]))
                atom4 = list(map(float, coords[int(atoms[3]) - 1]))
                dihedral_angle = misc.calculate_dihedral_angle(atom1,atom2, 
                                                               atom3,atom4)
                self.dihedrals_to_freeze.append((int(atoms[0]),int(atoms[1]), 
                                                 int(atoms[2]),int(atoms[3]), 
                                                 f"{dihedral_angle:.2f}"))
            excluded_numbers += [item for sublist in self.dihedrals_to_freeze \
                                 for item in sublist[:4]]

        excluded_numbers = list(set(excluded_numbers))
        self.excluded_numbers_str = ', '.join(map(str, excluded_numbers))
            
        bonds_text = set(item for item in excluded_numbers)

        included_numbers = [int(i) for i in range(1, startoutput.num_atom+1) \
                            if i not in bonds_text]
        self.included_numbers_str = misc.format_numbers_as_ranges(
                                    included_numbers)

        # Force constant
        self.force = re.search(r'^force constant(.+)', self.file_content,
                               flags=re.MULTILINE)
        if self.force:
            self.force = self.force.group(1).strip()
            force_value = float(self.force)

    def getactiveatoms(self):
        """Read active atoms from parameters file"""
        actats = re.search(r'^active atoms(.+)' , self.file_content,
                           flags=re.MULTILINE)
        if actats:
            self.activ_ats = [int(elt) for elt in actats.group(1).split()]
        else:
            print('Active atoms not specified in parameters file')
            sys.exit(1)

    def getrmsd(self):
        """Read RMSD thresholds"""
        rmsd_opt_def = 0.1
        rmsd_crest_def = 0.1
        #RMSD threshold
        RMSD = re.search(r'^rmsd threshold(.+)' , self.file_content,
                         flags=re.MULTILINE)
        if RMSD:
            if len(RMSD.group(1).split()) == 1:
                self.RMSD_crest = float(RMSD.group(1))
                self.RMSD_opt = rmsd_opt_def
            elif len(RMSD.group(1).split()) == 2: 
                self.RMSD_crest = float(RMSD.group(1).split()[0])
                self.RMSD_opt = float(RMSD.group(1).split()[1])
        else:
            self.RMSD_crest = rmsd_crest_def
            self.RMSD_opt = rmsd_opt_def

            
    def getopt(self):
        """Read informations about optimisation calculations"""
        software = re.search(r'^software (.+)', self.file_content,
                             flags=re.MULTILINE)
        if software:
            self.software = software.group(1).strip()
        else:
            self.software = "same"

        method = re.search(r'^method (.+)', self.file_content,
                           flags=re.MULTILINE)
        self.method = method.group(1).strip()

        basis_set = re.search(r'^basis set (.+)', self.file_content,
                              flags=re.MULTILINE)
        self.basis_set = basis_set.group(1).strip()
        pattern = re.compile('[^a-zA-Z0-9_]+')
        self.basis_name = pattern.sub('', self.basis_set) 

        dispersion = re.search(r'^dispersion (.+)', self.file_content,
                               flags=re.MULTILINE)
        if not dispersion:
            self.dispersion = ''
        else:
            self.dispersion = dispersion.group(1).strip()
            if self.dispersion == 'none': 
                self.dispersion = ''

        solvent = re.search(r'^opt solvent (.+)', self.file_content,
                            flags=re.MULTILINE)
        if not solvent:
            self.solvent = 'none'
        else:
            self.solvent = solvent.group(1).strip()

        cpus = re.search(r'^cpus (.+)', self.file_content,
                         flags=re.MULTILINE)
        if not cpus:
            self.cpus = 12
        else:
            self.cpus = cpus.group(1).strip()

        time = re.search(r'^time limit (.+)', self.file_content,
                         flags=re.MULTILINE)
        if not time:
            self.optime = '48:00:00'
        else:
            self.optime = time.group(1).strip()

        node = re.search(r'^nodes excluded (.+)', self.file_content,
                         flags=re.MULTILINE)
        if not node:
            self.node = 'none'
        else:
            self.node = node.group(1).strip()

        guidopt = re.search(r'^guided opt', self.file_content,
                            flags=re.MULTILINE)
        if guidopt:
            self.guidopt = True
        else:
            self.guidopt = False

        tsc = re.search(r'^TS screening (.+)' , self.file_content,
                        flags=re.MULTILINE)
        if tsc:
            self.tsc = tsc.group(1).strip()
        else:
            self.tsc = 'scalprod'

    def getadcalc(self):
        adcalc = re.search(r'^additional calculation (.+)', self.file_content,
                         flags=re.MULTILINE)
        if not adcalc:
            self.adcalc = "none"
        else:
            self.adcalc = [ac.strip() for ac in adcalc.group(1).split()]

        scope = re.search(r'^scope (.+)', self.file_content,
                          flags=re.MULTILINE)
        if not scope:
            self.scope = ["lowest"]
        else:
            self.scope = [sc.strip() for sc in scope.group(1).split()]

        if self.adcalc != "none":
            ad_method = re.search(r'^additional method (.+)', self.file_content,
                                  flags=re.MULTILINE)
            if not ad_method:
                self.ad_method = [self.method]
            else:       
                self.ad_method = [ad.strip() for ad in ad_method.group(1)\
                                                                .split()]

            ad_basis_set = re.search(r'^additional basis set (.+)', 
                                     self.file_content, flags=re.MULTILINE)
            if not ad_basis_set:
                self.ad_basis_set = [self.basis_set]
                self.ad_basis_name = [self.basis_name]
            else:       
                self.ad_basis_set = ad_basis_set.group(1).strip().split()
                pattern = re.compile('[^a-zA-Z0-9_]+')
                self.ad_basis_name = [pattern.sub('', bs) for bs in \
                                      self.ad_basis_set]

            ad_solvent = re.search(r'^additional solvent (.+)', 
                                   self.file_content, flags=re.MULTILINE)
            if not ad_solvent:
                self.ad_solvent = [self.solvent]
            else:       
                self.ad_solvent = ad_solvent.group(1).strip().split()

            ad_dispersion = re.search(r'^additional dispersion (.+)', 
                                   self.file_content, flags=re.MULTILINE)
            if not ad_dispersion:
                self.ad_dispersion = [self.dispersion]
            else:       
                self.ad_dispersion = ad_dispersion.group(1).strip().split()

            ad_software = re.search(r'^additional software(.+)', 
                                   self.file_content, flags=re.MULTILINE)
            if not ad_software:
                self.ad_software = [self.software]
            else:       
                self.ad_software = ad_software.group(1).strip().split()
        
    def getwoc(self):
        self.woc = re.search(r'^woc (.+)', self.file_content,
                             flags=re.MULTILINE)
        if not self.woc:
            self.woc = False
        else:
            if self.woc.group(1).strip() == 'True':
                self.woc = True
            else: 
                self.woc = False


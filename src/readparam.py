import sys
import re
import misc
import math

class Parameters:
    """Class that is used to read and store the informations in the parameters
    text file"""

    def __init__(self, file_name, startoutput=None):
        """Open the file and read """
        self.fn = file_name
        self.check_st = "unknown"
        with open(self.fn,'r') as file_brut:
            self.file_content = file_brut.read()
        with open(self.fn,'r') as file_brut:
            self.file_lines = file_brut.readlines()
        self.getcheck()
        self.getversion()
        self.getsolvent()
        self.getewin()
        self.getNCI()
        self.getrmsd()
        if startoutput != None:
            if startoutput.ts:
                self.getconstraints(startoutput)
        if self.fn.split('.')[-1] == 'txt':
            self.checkdft()
        elif self.fn.split('.')[-1] == 'tmp':
            self.getdft()
            self.getgeneral()
            if self.ref_freq:
                self.getactiveatoms()
        else:
            print(f"File {self.fn} should not be named like that")
            sys.exit(1)

    def getcheck(self):
        self.woc = re.search(r'WoC (.+)', self.file_content) 
        if not self.woc:
            self.woc = False
        else:
            if self.woc.group(1).strip() == 'True':
                self.woc = True
            else: 
                self.woc = False
        self.check = re.search(r'Check (.+)', self.file_content) 
        if not self.check:
            self.check = 'none'
        elif self.check == 'none':
            self.check = self.check.group(1).strip() 
        else:
            try:
                self.check = int(self.check.group(1).strip())
            except:
                print('Check should be none or an integer')
                sys.exit(1)
        
    def getversion(self):
        """Read which version of CREST to use"""
        self.version_crest = re.search(r'CREST version(.+)', self.file_content) 
        if self.version_crest:
            self.version_crest = self.version_crest.group(1).strip() 
            if self.version_crest == 'default': 
                self.version_crest = 'crest'
            elif self.version_crest == '3.0':
                self.version_crest = 'crest3'
            elif self.version_crest == 'continous release':
                self.version_crest = 'crest_continous_release'
            else:
                print('Please provide a suitable CREST version.')
                sys.exit()
        else:
            print('Unable to find CREST version information in the parameters'
                  ' file.')
            sys.exit()

    def getewin(self):
        """Read the size of the energy window, if none, take the default value
        6 kcal/mol"""
        self.ewin_crest = re.search(r'EWIN(.+)', self.file_content) 
        
        if self.ewin_crest:
            self.ewin_crest = float(self.ewin_crest.group(1).strip())
        else:
            self.ewin_crest = 6
    def getsolvent(self):
        self.solvent_crest = re.search(r'CREST solvent(.+)', self.file_content)

        if self.solvent_crest:
            self.solvent_crest = self.solvent_crest.group(1).strip()

            if self.solvent_crest.lower() == 'none' :
                self.solvent_crest = ''
            elif re.search(r'--alpb', self.solvent_crest):
                pass
            else:
                print('Please provide a suitable CREST solvent.')
                sys.exit()
        else:
            print('Unable to find CREST solvent information in the parameters '
                  'file.')
            sys.exit()

    def getNCI(self):
        self.nci = re.search(r'NCI(.+)', self.file_content)

        if self.nci:
            self.nci = self.nci.group(1).strip() 

            if self.nci.lower() == 'none':
                self.nci = ''
            elif re.search(r'selected', self.nci):
                if re.search(r'selected with scale factor', self.nci):
                    nci_factor = re.search(r'selected with scale factor(.+)', 
                                           self.nci).group(1).strip()
                    try:
                        self.nci_value = float(nci_factor)
                        self.nci = f'--nci --wscal {self.nci_value}'
                    except ValueError:
                        print('Please provide a suitable NCI scale factor.')
                        sys.exit()  
                else:
                    self.nci = '--nci'
            else:
                print('Please provide a suitable CREST NCI.')
                sys.exit()
        
        else:
            print('Unable to find CREST NCI information in the parameters '
                  'file.')
            sys.exit() 

    def getconstraints(self,startoutput):
        bonds = re.search('Bond constrained(.+)', self.file_content)
        angles = re.search('Angle constrained(.+)', self.file_content)
        dihedrals = re.search('Dihedral angle constrained (.+)', 
                              self.file_content)
        
        if bonds or angles or dihedrals: 
            self.bonds_to_freeze = []
            self.angles_to_freeze = []
            self.dihedrals_to_freeze = []
            lines = startoutput.coordinates
        else:
            print('Unable to find constrained bond, angle or dihedral angle'
                  ' information in the parameters file.')
            sys.exit()
        
        if bonds:
            for line in self.file_lines:
                if 'Bond constrained' in line:
                    atoms = line.split()[2:]
                    if len(atoms) == 2:
                        atom1 = list(map(float, lines[int(atoms[0]) - 1]))
                        atom2 = list(map(float, lines[int(atoms[1]) - 1]))
                        distance = math.dist(atom1, atom2)
                        self.bonds_to_freeze.append((int(atoms[0]), 
                                                     int(atoms[1]), 
                                                    "{:.3f}".format(distance)))
                    else:
                        print('Please provide a suitable constrained bond.')
                        sys.exit()

        if angles:
            for line in parameters_lines:
                if 'Angle constrained' in line:
                    atoms = line.split()[2:]
                    if len(atoms) == 3:
                        atom1 = list(map(float, lines[int(atoms[0]) - 1]))
                        atom2 = list(map(float, lines[int(atoms[1]) - 1]))
                        atom3 = list(map(float, lines[int(atoms[2]) - 1]))
                        angle = misc.calculate_angle(atom1, atom2, atom3)
                        self.angles_to_freeze.append((int(atoms[0]), 
                                                      int(atoms[1]), 
                                                      int(atoms[2]), 
                                                 "{:.2f}".format(angle)))
                    else:
                        print('Please provide a suitable constrained angle.')
                        sys.exit()

        if dihedrals:
            for line in parameters_lines:
                if 'Dihedral angle constrained' in line:
                    atoms = line.split()[3:]
                    if len(atoms) == 4:
                        atom1 = list(map(float, lines[int(atoms[0]) - 1]))
                        atom2 = list(map(float, lines[int(atoms[1]) - 1]))
                        atom3 = list(map(float, lines[int(atoms[2]) - 1]))
                        atom4 = list(map(float, lines[int(atoms[3]) - 1]))
                        dihedral_angle = misc.calculate_dihedral_angle(atom1, 
                                                                       atom2, 
                                                                       atom3, 
                                                                       atom4)
                        self.dihedrals_to_freeze.append((int(atoms[0]), 
                                                         int(atoms[1]), 
                                                         int(atoms[2]), 
                                                         int(atoms[3]), 
                                             "{:.2f}".format(dihedral_angle)))
                    else:
                        print('Please provide a suitable constrained dihedral'
                              ' angle.')
                        sys.exit()

        excluded_numbers_bonds =  []
        excluded_numbers_bonds = [item for sublist in self.bonds_to_freeze for \
                                  item in sublist[:2]]

        excluded_numbers_angle = []
        excluded_numbers_angle = [item for sublist in self.angles_to_freeze \
                                  for item in sublist[:3]]

        excluded_numbers_dihedral = []
        excluded_numbers_dihedral = [item for sublist in \
                                    self.dihedrals_to_freeze for item in \
                                    sublist[:4]]

        excluded_numbers = []
        excluded_numbers = excluded_numbers_bonds + excluded_numbers_angle \
                         + excluded_numbers_dihedral
        excluded_numbers = list(set(excluded_numbers))
        self.excluded_numbers_str = ', '.join(map(str, excluded_numbers))
            
        bonds_text = set(item for item in excluded_numbers)

        included_numbers = [int(i) for i in range(1, startoutput.num_atom+1) \
                            if i not in bonds_text]
        self.included_numbers_str = misc.format_numbers_as_ranges(
                                    included_numbers)

        # Force constant
        self.force = re.search(r'Force constant(.+)', self.file_content)

        if self.force:
            self.force = self.force.group(1).strip()

            try:
                force_value = float(self.force)
            except ValueError:
                print('Please provide a suitable force constant.')
                sys.exit()
        
        else:
            print('Unable to find force constant information in the parameters'
                  ' file.')
            sys.exit()      

    def getactiveatoms(self):
        actats = re.search(r'Active atoms(.+)' , self.file_content)
        if actats:
            self.activ_ats = [int(elt) for elt in actats.group(1).split()]
        else:
            print('Active atoms not specified in parameters file')
            sys.exit(1)

    def getrmsd(self):
        #RMSD threshold
        RMSD = re.search(r'RMSD threshold(.+)' , self.file_content)
    
        if RMSD:
            if len(RMSD.group(1).split()) == 1:
                RMSD_crest = RMSD.group(1)
                RMSD_dft = RMSD.group(1)
                try:
                    self.RMSD_crest = float(RMSD_crest)
                    self.RMSD_dft = float(RMSD_dft)
                except ValueError:
                    print('Please provide a suitable RMSD threshold.')
                    sys.exit(1)
            elif len(RMSD.group(1).split()) == 2: 
                RMSD_crest = RMSD.group(1).split()[0]
                RMSD_dft = RMSD.group(1).split()[1]
                try:
                    self.RMSD_crest = float(RMSD_crest)
                except ValueError:
                    print('Please provide a suitable RMSD threshold.')
                    sys.exit(1)
                try:
                    self.RMSD_dft = float(RMSD_dft)
                except ValueError:
                    print('Please provide a suitable RMSD threshold.')
                    sys.exit(1)
            else:
                print('Too many RMSD thresholds in parameter file')
                sys.exit(1) 
        
        else:
            print('Unable to find RMSD threshold information in the parameters'
                  ' file.')
            sys.exit() 

    def checkdft(self):
        compsofts = ['gaussian','orca']
        addcalcs = ['irc', 'nbo', 'cdft', 'none']
        #software
        self.software = re.search(r'Software(.+)', self.file_content)
        if self.software:
            if self.software.group(1).strip().lower() not in compsofts:
                print(f'Software is unknown, it must be one of {compsofts}.')
                sys.exit()   
        else:
            self.software = "same"
                
        #functionnal
        functional = re.search(r'Functional(.+)', self.file_content)
    
        if functional:
            pass
        else:
            print('Unable to find functional information in the parameters '
                  'file.')
            sys.exit()   
    
        #base
        base_1 =  re.search(r'First base(.+)', self.file_content)
        base_2 = re.search(r'Unique base(.+)', self.file_content)
    
        if base_1 or base_2:
            pass
        else:
            print('Unable to find basis set information in the parameters '
                  'file.')
            sys.exit()
    
        #dispersion
        dispersion = re.search(r'Dispersion(.+)', self.file_content)
    
        if dispersion:
            pass
        else:
            print('Unable to find dispersion information in the parameters '
                  'file.')
            sys.exit() 
    
        # solvent
        solvent = re.search(r'DFT solvent(.+)', self.file_content)
        if solvent:
            pass
        else:
            print('Unable to find solvent information in the parameters file.')
            sys.exit()
    
        #additional calculation
        self.choice = re.search(r'Additional calculation (.+)', 
                                self.file_content)
        self.choice = self.choice.group(1).strip().lower()
        if not self.choice:
            print("Unable to find additional calculation information in the "
                  "parameters file.")
            sys.exit(1) 
        elif self.choice not in addcalcs:
            print("Additional calculation should be one of the follow:\n")
            for ac in addcalcs:
                print(f"- {ac}")
            sys.exit(1) 
            
        #node for calculation
        node = re.search(r'Excluding nodes (.+)', self.file_content)
        if node:
            pass
        else:
            print("Unable to find node excluded information in the parameters "
                   "file.")
            sys.exit() 
            
    def getdft(self):
        software = re.search(r'Software (.+)', self.file_content)
        if software:
            self.software = software.group(1).strip()
        else:
            self.software = "same"

        functional = re.search(r'Functional (.+)', self.file_content)
        self.functional = functional.group(1).strip()

        base_1 =  re.search(r'First base (.+)', self.file_content)
        if base_1:
            self.base_1 = base_1.group(1).strip()
            base_2 = re.search(r'Second base (.+)', self.file_content)
            self.base_2 = base_2.group(1).strip()
        else:
            base_1 = re.search(r'Unique base (.+)', self.file_content)
            self.base_1 = base_1.group(1).strip()
            self.base_2 = self.base_1

        #pattern = re.compile('[\W_]+')
        pattern = re.compile('[^a-zA-Z0-9_]+')
        self.base_name_1 = pattern.sub('', self.base_1)
        self.base_name_2 = pattern.sub('', self.base_2)

        dispersion = re.search(r'Dispersion (.+)', self.file_content)
        self.dispersion = dispersion.group(1).strip()
        if self.dispersion.lower() == 'none': 
            self.dispersion = ''

        solvent = re.search(r'DFT solvent (.+)', self.file_content)
        self.solvent = solvent.group(1).strip()
        if self.solvent.lower() == 'none':
            self.solvent = ''

        self.choice = re.search(r'Additional calculation (.+)', 
                                self.file_content)
        self.choice = self.choice.group(1).strip().lower()
        if self.choice == 'none' :
            self.choice = ''

        node = re.search(r'Excluding nodes (.+)', self.file_content)
        self.node = node.group(1)
        if self.node.lower() == 'none':
            self.node = ''
        
    def getgeneral(self):
        experience_number = re.search(r'Experience name (.+)',
                                           self.file_content)
        self.experience_number = experience_number.group(1).strip()
    
        charge = re.search(r'Charge (-?\d+)', self.file_content)
        self.charge = charge.group(1).strip()

        multiplicity = re.search(r'Multiplicity (-?\d+)', self.file_content)
        self.multiplicity = multiplicity.group(1).strip()

        ref_freq = re.search(r'First frequency (.+)', self.file_content)
        if ref_freq.group(1).strip().lower() == 'none':
            self.ref_freq = ''
        else:
            self.ref_freq = int(ref_freq.group(1))

        check_st = re.search(r'Loop number (-?\d+)', self.file_content)
        self.check_st = int(check_st.group(1).strip())


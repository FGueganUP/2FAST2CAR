import sys
import cclib
import re
import os
fastcar_dir = "/home/hgeindre/ProjectFASTCAR/FASTCAR2"
sys.path.append(fastcar_dir)
import filemanager as fm

class Output(cclib.parser.data.ccData):
    """Class that is used to read and store information from the outputfile of 
    a quantum chemistry software used as starting point."""

    def __init__(self, filename):
        """Open the file, identify software, check termination and read 
        charge, multiplicity and geometry"""
        conv = 0.036749308136649
        # Read the file with cclib (https://cclib.github.io/index.html)
        parser = cclib.io.ccopen(filename)
        data = parser.parse()
        self.file_name = filename
        # Name of the file stripped of extension and path
        self.base_name = os.path.splitext(self.file_name.split('/')[-1])[0]
        #print(data.metadata)
        if '_lowest' in self.base_name:
            self.base_name = self.base_name.strip('_lowest')
        self.software = data.metadata['package']
        #self.normterm = data.optdone
        self.normterm = data.metadata['success']
        if 'wall_time' in data.metadata.keys():
            self.realtime = data.metadata['wall_time']
        else:
            self.realtime = 'none'
        if 'cpu_time' in data.metadata.keys():
            self.cputime = data.metadata['cpu_time']
        else:
            self.cputime = 'none'
        if not self.normterm:
            return
        if data.metadata['methods'][0] == 'DFT':
            if 'functional' in data.metadata.keys():
                self.method = data.metadata['functional'].lower()
            else:
                self.method = 'dft'
        else:
            self.method = data.metadata['methods'][0].lower()
        self.basis_set = data.metadata['basis_set'].lower()
        if 'solvent_name' in data.metadata:
            self.solvent = data.metadata['solvent_name'].lower()
        else: 
            self.solvent = 'none'
        self.charge = data.charge
        self.mult = data.mult
        self.last_SCF = data.scfenergies[-1]*conv
        self.energies = [nrj*conv for nrj in data.scfenergies]
        self.num_atom = data.natom
        self.coordinates = data.atomcoords[-1]
        self.atnums = data.atomnos
        if hasattr(data, 'vibfreqs'):
            self.freqs = data.vibfreqs
            if self.freqs[0] < 0:
                self.ts = True
            else:
                self.ts = False
            self.normvec = data.vibdisps
        if hasattr(data, 'atomcharges'):
            self.charges = data.atomcharges

    #def idsoft(self):
    #    """Identify the quantum chemistry software used in the outputfile.
    #    Compatible softwares: Gaussian, Orca"""
    #    if "Entering Gaussian System" in self.file_content:
    #        self.software = "gaussian"
    #    elif "O   R   C   A" in self.file_content:
    #        self.software = "Orca"
    #    else: 
    #        print("Unknown software")
    #        sys.exit(1)

    #def checkterm(self):
    #    """Check if the calculation ended properly.
    #    Compatible softwares : Gaussian, Orca"""
    #    if self.software == "gaussian" and\
    #       "Normal termination of Gaussian 16" in self.file_content:
    #        self.normterm = True
    #    elif self.software == "Orca" and\
    #         "ORCA TERMINATED NORMALLY" in self.file_content:
    #        self.normterm = True
    #    else:
    #        self.normterm = False
    #
    #def read_chamul(self):
    #    """Read charge and multiplicity in the outputfile.
    #    Compatible softwares: Gaussian, Orca"""
    #    if self.software == "gaussian":
    #        charge_match = re.search(r'Charge\s*=\s*(-?\d+)\s*', 
    #                                 self.file_content)
    #        self.charge = int(charge_match.group(1))

    #        multiplicity_match = re.search(r'Multiplicity\s*=\s*(-?\d+)\s*', 
    #                                       self.file_content)
    #        self.multiplicity = int(multiplicity_match.group(1))
    #    elif self.software == "Orca":
    #        charge_match = re.search(r'Total Charge\s*Charge.*(-?\d+)\s*', 
    #                                 self.file_content)
    #        self.charge = int(charge_match.group(1))

    #        multiplicity_match = re.search(r'Multiplicity\s*Mult.*(-?\d+)\s*', 
    #                                       self.file_content)
    #        self.multiplicity = int(multiplicity_match.group(1))

    #def readtheorylvl(self):
    #    """Read calculation method and basis set
    #    Compatible software: Gaussian"""
    #    if self.software == "gaussian":
    #        self.method = re.findall("SCF Done:  E[(](.*)[)]", 
    #                                 self.file_content)
    #        self.method = self.method[0].lower()
    #        self.basis_set = re.findall("Standard basis: (.*?)\n", 
    #                                    self.file_content)
    #        self.basis_set = self.basis_set[0].split()[0].lower()
    #        
    #def readsolvent(self):
    #    """Read solvent used in calculation
    #    Compatible software: Gaussian"""
    #    if self.software == "gaussian":
    #        self.solvent = re.search("scrf=[(](.*?)[)]", 
    #                                 self.file_content.lower())
    #        if self.solvent:
    #            self.solvent = self.solvent.group()
    #        else:
    #            self.solvent = ''

    #def readnrj(self):
    #    """Read the last SCF energy. Compatible software : Gaussian"""
    #    if self.software == "gaussian":
    #        self.last_SCF = re.findall("SCF Done:  E[(].*", self.file_content)
    #        self.last_SCF = float(self.last_SCF[-1].strip().split()[4]) \
    #            if self.last_SCF else None

    #def readallnrjs(self):
    #        """Read all energies, with ZPE, thermal energy, thermal enthalpy,
    #        free energy. Compatible software : Gaussian"""
    #        ElecE = re.findall("Sum of electronic and zero-point Energies=.*",
    #                           self.file_content)
    #        ElecE = ElecE[0].strip().split()[6] if ElecE else None
    #        ThermalE = re.findall("Sum of electronic and thermal Energies=.*",
    #                              self.file_content)
    #        ThermalE = ThermalE[0].strip().split()[6] if ThermalE else None
    #        
    #        EnthalpiesE = re.findall("Sum of electronic and thermal Enthalpies"
    #                                 "= .*", self.file_content)
    #        EnthalpiesE = EnthalpiesE[0].strip().split()[6] if EnthalpiesE \
    #                      else None
    #        
    #        FreeE = re.findall("Sum of electronic and thermal Free Energies=.*"
    #                           ,self.file_content)
    #        FreeE = FreeE[0].strip().split()[7] if FreeE else None

    #        self.allnrjs = ((self.file_name,self.last_SCF,ElecE,ThermalE,
    #                         EnthalpiesE,FreeE,' '))
    #    
    #def extractgeo(self):
    #    """Extract last coordinates found in file.
    #    Compatible softwares : Gaussian, Orca"""
    #    ecs = fm.readepc(keys='letter')
    #    if self.software == "gaussian":
    #        pattern = r'Standard orientation:' + '(.*?)' + 'Rotational constants'
    #        match = re.findall(pattern, self.file_content, re.DOTALL)
    #        self.coordinates = list()
    #        self.atnums = list()
    #        for line in match[-1].split('\n')[5:-2]:
    #            self.atnums.append(int(line.split()[1]))
    #            self.coordinates.append([float(elt) for elt in line.split()[3:]])
    #        self.num_atom = len(self.coordinates)
    #    elif self.software == "Orca":
    #        # !!!!!!!!!!!!!!! A REFAIRE !!!!!!!!!!!!!!!!!  
    #        pattern = r'CARTESIAN COORDINATES \(ANGSTROEM\)' + '(.*?)'\
    #                + r'CARTESIAN COORDINATES \(A.U.\)'
    #        match = re.finditer(pattern, self.file_content, re.DOTALL)
    #        for m in match:
    #            self.coordinates = re.sub(r'-{2,}','',m.group(1))
    #            self.coordinates = re.sub(r'(?m)^\s$\n?','',self.coordinates)
    #
    #        self.num_atom = len(self.coordinates.split('\n')) - 1

    #def extractgeo_old(self):
    #    """Extract last coordinates found in file.
    #    Compatible softwares : Gaussian, Orca"""
    #    ecs = fm.readepc(keys='letter')
    #    if self.software == "gaussian":
    #        pattern = r'\\' + '(.*)' + r'\\'
    #        match = re.search(pattern, self.file_content, re.DOTALL)
    #        coordinates = match.group(1)
    #        coordinates = re.sub(r'[\r\n\s]+', '', coordinates)
    #        pattern = r'\\' + f'{self.charge},{self.multiplicity}'\
    #                  + r'\\(.*?)\\Version'
    #        match = re.search(pattern, coordinates)
    #        coordinates = match.group(1)
    #        self.coordinates = list()
    #        self.atnums = list()
    #        for line in (coordinates.split('\\')):
    #            if not line:
    #                continue
    #            self.coordinates.append([float(elt) for elt in line.split(',')[1:]])
    #            self.atnums.append(ecs[line.split(',')[0]][0])

    #        self.num_atom = len(self.coordinates)
    #    elif self.software == "Orca":
    #        # !!!!!!!!!!!!!!! A REFAIRE !!!!!!!!!!!!!!!!!  
    #        pattern = r'CARTESIAN COORDINATES \(ANGSTROEM\)' + '(.*?)'\
    #                + r'CARTESIAN COORDINATES \(A.U.\)'
    #        match = re.finditer(pattern, self.file_content, re.DOTALL)
    #        for m in match:
    #            self.coordinates = re.sub(r'-{2,}','',m.group(1))
    #            self.coordinates = re.sub(r'(?m)^\s$\n?','',self.coordinates)
    #
    #        self.num_atom = len(self.coordinates.split('\n')) - 1

    #def extractcoord(self):
    #    """Extract redundant internal coordinates for TS.
    #    Compatible software : Gaussian"""
    #    start_phrase = 'Redundant internal coordinates found in file.  '\
    #                   '(old form).'
    #    end_phrase = 'Recover connectivity data from disk.'
    #    pattern = re.escape(start_phrase) + '(.*?)' + re.escape(end_phrase)
    #    match = re.search(pattern, self.file_content, re.DOTALL)
    #    coordinates_text = match.group(1)
    #    coordinates_text = coordinates_text.replace(",0,", "         ")
    #    self.coord_text = coordinates_text.replace(",", "        ")  
    #    self.coord_text = self.coord_text.lstrip()

    #def readfreqs(self):
    #    """Read and store in a list the 3n-6 first (Gaussian) or last (Orca) 
    #    frequencies found in the output. Software compatibles : Gaussian, 
    #    Orca"""
    #    if self.software == "gaussian":
    #        pattern = r'\s*Frequencies --\s*(.*)'
    #        self.freqs = list()
    #        if not re.search(pattern, self.file_content):
    #            self.freqs.append(1)
    #        else:
    #            match = re.finditer(pattern, self.file_content)
    #            if self.num_atom > 2:
    #                nfreqs = 3*self.num_atom - 6
    #            else:
    #                nfreqs = 3*self.num_atom - 5
    #            for ms in match:
    #                for freq in ms.group(1).split():
    #                    self.freqs.append(float(freq))
    #                if len(self.freqs) == nfreqs:
    #                    break
    #    elif self.software == "Orca":
    #        pattern = r'\s*VIBRATIONAL FREQUENCIES\s*\n\s*-+\n(.*?)\n-+\n'
    #        match = re.finditer(pattern, self.file_content, re.DOTALL)
    #        if match == None:
    #            print("No frequencies calculation in the outputfile")
    #            sys.exit(1)
    #        self.freqs = list()
    #        if self.num_atom > 2:
    #            nfreqs = 3*self.num_atom - 6
    #        else:
    #            nfreqs = 3*self.num_atom - 5
    #        for ms in reversed(list(match)):
    #            for line in ms.group(1).split('\n'):
    #                if not line:
    #                    continue
    #                freq = float(line.split()[1])
    #                if freq != 0.0:
    #                    self.freqs.append(freq)
    #            if len(self.freqs) == nfreqs:
    #                break
    #    
    #    if self.freqs[0] < 0:
    #        self.ts = True
    #    else:
    #        self.ts = False

    #def readnormvec(self):
    #    """Read vectors of normal modes
    #    Compatible softwares: Gaussian"""
    #    pattern = r'\s*Frequencies --\s*(.*)'
    #    match = re.search(pattern, self.file_content)
    #    self.normvec = list()
    #    for line in self.file_content[match.end():]\
    #                .split('\n')[5:5+self.num_atom]:
    #        sline = line.split()
    #        self.normvec.append((int(sline[1]),
    #                            [float(sline[i]) for i in range(2,5)]))

    def readcharges(self):
        """Read natural population analysis. 
        Compatible software : Gaussian"""
        self.charges = {}

        self.charges['npa'] = []
        npa_start_phrase = ' Summary of Natural Population Analysis:'
        npa_end_phrase = ' ==================================================='
        npa_pattern = re.escape(npa_start_phrase) + '(.*?)' \
                    + re.escape(npa_end_phrase)
        npa_match = re.search(npa_pattern, self.file_content, re.DOTALL)
        npa_text = npa_match.group(1)
        lines = npa_text.split('\n')
        lines = lines[6:-1]
        for line in lines:
            self.charges['npa'].append((line.split()[0],
                                        float(line.split()[2])))

        self.charges['mul'] = []
        mul_start_phrase = " Mulliken charges"
        mul_end_phrase = " Sum of"
        mul_pattern = mul_start_phrase + '(.*?)' + mul_end_phrase
        mul_match = re.search(mul_pattern, self.file_content, re.DOTALL)
        mul_text = mul_match.group(1)
        lines = mul_text.split('\n')
        lines = lines[2:-1]
        for line in lines:
            self.charges['mul'].append((line.split()[1],
                                        float(line.split()[2])))
    
        #self.charges['apt'] = []
        #apt_start_phrase = " APT charges"
        #apt_end_phrase = " Sum of"
        #apt_pattern = apt_start_phrase + '(.*?)' + apt_end_phrase
        #apt_match = re.search(apt_pattern, self.file_content, re.DOTALL)
        #apt_text = apt_match.group(1)
        #lines = apt_text.split('\n')
        #lines = lines[2:-1]
        #for line in lines:
        #    self.charges['apt'].append((line.split()[1],
        #                                float(line.split()[2])))

        self.charges['hir'] = []
        hir_start_phrase = " Hirshfeld charges"
        hir_end_phrase = " Tot"
        hir_pattern = hir_start_phrase + '(.*?)' + hir_end_phrase
        hir_match = re.search(hir_pattern, self.file_content, re.DOTALL)
        hir_text = hir_match.group(1)
        lines = hir_text.split('\n')
        lines = lines[2:-1]
        for line in lines:
            self.charges['hir'].append((line.split()[1],
                                        float(line.split()[2])))
    
    def readnbo(self):
        """Read NBO. Compatible software : Gaussian"""
        self.NBO = []
        start_phrase = 'NATURAL BOND ORBITALS (Summary):'
        end_phrase = 'NBO analysis completed in'
        pattern = re.escape(start_phrase) + '(.*?)' + re.escape(end_phrase)
        match = re.search(pattern, self.file_content, re.DOTALL)
        coordinates_text = match.group(1)
        lines_with_BD = [line for line in coordinates_text.strip().split('\n')\
                          if re.search(r'BD\*', line)]
        for item in lines_with_BD:
            item = item.split()
            bond_info = (item[3],item[4],item[5],item[6])
            bond_info = ''.join(bond_info)
            bond_value = item[8]
            self.NBO.append((bond_info,bond_value))
        self.NBO.sort(key=lambda item: item[1])


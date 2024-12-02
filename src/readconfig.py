import sys
import os

class Cluster:
    """Class that is used to read and store the informations about the 
    cluster"""

    def __init__(self,fastcar_dir):
        """Open the file and read """
        self.fn = fastcar_dir + "/config.txt"
        if not os.path.isfile(self.fn):
            print(f'File {self.fn} missing.')
            sys.exit()
        with open(self.fn,'r') as file_brut:
            self.file_lines = file_brut.readlines()
        for line in self.file_lines:
            if line.startswith('#'):
                continue
            elif line.startswith("g16_path"):
                self.g16_path = line.split("=")[1].strip()
            elif line.startswith("g16C_path"):
                self.g16C_path = line.split("=")[1].strip()
        if not hasattr(self, "g16_path"):
            print(f'Variable g16_path missing in {self.fn}')
            sys.exit()
    
class Keywords:
    """Class that is used to read and store the informations about the 
    keywords use in the quantum chemistry softwares"""
    def __init__(self,fastcar_dir):
        """Open the file, and store informations in a dict """
        self.keywords = dict()
        self.fn = fastcar_dir + "/keywords.txt"
        if not os.path.isfile(self.fn):
            print(f'File {self.fn} missing.')
            sys.exit()
        with open(self.fn,'r') as file_brut:
            self.file_lines = file_brut.readlines()
        for line in self.file_lines:
            if line.startswith('#'):
                continue
            elif line:
                soft = line.split(';')[0].strip().lower()
                calctype = line.split(';')[1].strip().lower()
                kwds = line.split(';')[2].strip("\n")
                if not soft in self.keywords:
                    self.keywords[soft] = dict()
                self.keywords[soft][calctype] = kwds.lower()


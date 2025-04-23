import sys
import os

class Cluster:
    """Class that is used to read and store the informations about the 
    cluster"""

    def __init__(self,fastcar_dir):
        """Open the file and read """
        self.fn = fastcar_dir + "/config.txt"
        self.default = dict()
        self.path = dict()
        if not os.path.isfile(self.fn):
            print(f'File {self.fn} missing.')
            sys.exit()
        with open(self.fn,'r') as file_brut:
            self.file_lines = file_brut.readlines()
        for line in self.file_lines:
            if line.startswith('#'):
                continue
            elif line.startswith("default"):
                soft = line.split('=')[0].split('_')[1].strip()
                self.default[soft] = line.split("=")[1].strip()
            elif line.startswith("path"):
                soft = line.split('=')[0].split('_')[1].strip()
                self.path[soft] = line.split("=")[1].strip()
            elif line.startswith("sub_command"):
                self.subco = line.split("=")[1].strip()
    
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


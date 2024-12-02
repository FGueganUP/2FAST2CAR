#!/usr/bin/env python3
import sys
import os
import glob
import re

fastcar_dir = os.path.dirname(os.path.realpath(__file__))

pattern = re.compile(r'^fastcar_dir = (.*)\n', re.MULTILINE)
for pyfile in (glob.glob('*.py')):
    mod = False
    if pyfile == os.path.basename(__file__):
        continue
    with open(pyfile,'r') as f:
        file_content = f.read()
        match = pattern.search(file_content)
        if match:
            mod = True
            new_content = file_content.replace(match.group(1),
                                               f'"{fastcar_dir}"')
    if mod:
        with open(pyfile, 'w') as f:
            f.write(new_content)
        
#Add the fastcar_dir in path

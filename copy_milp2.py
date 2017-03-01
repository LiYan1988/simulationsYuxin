# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 11:43:00 2017

@author: yx4vf

This file copies milp_before_cluster\milp2.py to everywhere there is milp2.py
after changes are made to milp_before_cluster\milp2.py
"""

import os
import shutil

for dirpath, dirnames, filenames in os.walk(os.path.curdir):
    if dirpath!='.\milp_before_cluster' and 'milp2.py' in filenames:
        print(dirpath)        
        shutil.copy('.\milp_before_cluster\milp2.py', dirpath)
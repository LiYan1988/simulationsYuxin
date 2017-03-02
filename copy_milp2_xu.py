# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 10:31:06 2017

@author: yx4vf

This file copies milp_before_cluster\milp2_xu.py to everywhere there is 
milp2_xu.py after changes are made to milp_before_cluster\milp2_xu.py
"""



import os
import shutil

for dirpath, dirnames, filenames in os.walk(os.path.curdir):
    if dirpath!='.\milp_before_cluster' and 'milp2_xu.py' in filenames:
        print(dirpath)        
        shutil.copy('.\milp_before_cluster\milp2_xu.py', dirpath)
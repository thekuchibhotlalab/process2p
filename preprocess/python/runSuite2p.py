# -*- coding: utf-8 -*-
"""
Created on Fri Aug 21 13:37:35 2020

@author: zzhu34
"""

import sys
import os  
sys.path.insert(0, 'C:/Users/zzhu34/Documents/gitRep/suite2p')
import suite2p
from suite2p import run_s2p



ops = run_s2p.default_ops() # populates ops with the default options
ops['nplanes'] = 2
ops['nchannels'] = 1 # Only for green channel animal
ops['fs'] = 15.63
ops['tau'] = 0.7
ops['functional_chan'] = 1
ops['save_mat'] = True
ops['align_by_chan'] : 1
ops['look_one_level_down'] = 1
ops['diameter'] = 10
print(ops)
days = [5,6];

data_folder = 'D:/labData/zz018h5' + '/' + 'zz018_016_suite2p'

db = {
        'h5py': 'D:/labData/zz018h5' + '/' + 'zz018_016_000.h5', # a single h5 file 
        'h5py_key': ['data'], # list of keys to use (they will be extracted in the order you give them
        'data_path': data_folder
        }

#opsEnd=run_s2p.run_s2p(ops=ops,db=db)

#for day in days: 
#    data_folder = 'K:/Jenni/se063/day' + str(day)
#    listfiles = os.listdir(data_folder)
#    db = {
#        'h5py': 'K:/Jenni/se063/day' + str(day) + '/' + listfiles[0], # a single h5 file 
#        'h5py_key': ['data'], # list of keys to use (they will be extracted in the order you give them
#        'data_path': data_folder
#        }

    # run one experiment
    

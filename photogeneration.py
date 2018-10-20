# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19, 2018

@author: Tim
"""

import numpy as np

def get_photogeneration(params):
    

    try:
        gen_file = open(params.gen_rate_file_name, "r")
    except:
        print(f"Unable to open file{params.gen_rate_file_name}")
            
    photogen_rate = gen_file.read()
    
    photogen_rate = params.photogen_scaling * photogen_rate/np.max(photogen_rate)
    
    gen_file.close()
            
  
            
        
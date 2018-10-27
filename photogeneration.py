# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19, 2018

@author: Timofey  Golubev

This just contains the function for reading photogeneration rate from a generation rate data file.
"""

import numpy as np

def get_photogeneration(params):
    '''
    Reads photogeneration rate from an input file.
    Inputs: 
        params: Params object which contains several necessary parameters such as the name of generation
                rate file as well as the photogeneration scaling to use. The photogeneration scaling
                is determined by finding a value which will result in the correct short-circuit current.
    '''
    
    try:
        gen_file = open(params.gen_rate_file_name, "r")
    except:
        print(f"Unable to open file{params.gen_rate_file_name}")
            
    photogen_rate = np.loadtxt(params.gen_rate_file_name)  
    #photogen_rate = np.ones(params.num_cell+1)
    
    photogen_rate = params.Photogen_scaling * photogen_rate/np.max(photogen_rate)
    
    gen_file.close()
    
    return photogen_rate
            
  
            
        
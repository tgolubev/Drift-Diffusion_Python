# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19, 2018

@author: Tim
"""

import numpy as np
import constants as const
from numba import jit

class Recombo():
    
    def __init__(self, params):
        self.R_Langevin = np.zeros(params.num_cell)
    
    @jit  
    def compute_R_Langevin(self, R_Langevin, n, p, N, k_rec, n1, p1):
        
        R_Langevin[1:] = k_rec*(N*N*n[1:]*p[1:] - n1*p1)
        
        # negative recombination values are unphysical
        for val in R_Langevin:
            if val < 0: 
                val = 0
                
        return R_Langevin
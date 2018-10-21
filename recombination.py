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
        #self.k_rec = params.k_rec
        self.R_Langevin = np.zeros(params.num_cell)
        #self.E_trap = params.active_VB + params.E_gap/2.0  #trap assisted recombo is most effective when trap is located mid-gap--> take VB and add 1/2 of bandgap
        #self.n1 = params.N_LUMO*np.exp(-(params.active_CB - self.E_trap)/const.Vt)
        #self.p1 = params.N_HOMO*np.exp(-(self.E_trap - params.active_VB)/const.Vt)
    
    @jit    
    def compute_R_Langevin(self, R_Langevin, n, p, N, k_rec, n1, p1):
        
        #N = params.N
        #k_rec = params.k_rec
        #n1 = params.n1
        #p1 = params.p1
        
        for i in range(1, len(p)):
            R_Langevin[i] = k_rec*(N*N*n[i]*p[i] - n1*p1)
            
            if R_Langevin[i] < 0: 
                R_Langevin[i] = 0
                
        return R_Langevin
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19, 2018

@author: Tim
"""

import numpy as np
import constants as const

class Recombo():
    
    def __init__(self, params):
        self.k_rec = params.k_rec
        self.R_langevin = np.zeros(params.num_cell)
        self.E_trap = params.active_VB + params.E_gap/2.0  #trap assisted recombo is most effective when trap is located mid-gap--> take VB and add 1/2 of bandgap
        self.n1 = params.N_LUMO*np.exp(-(params.active_CB - self.E_trap)/const.Vt)
        self.p1 = params.N_HOMO*np.exp(-(self.E_trap - params.active_VB)/const.Vt)
        
    def compute_R_Langevin(self, params, n, p):
        
        for i in range(1, len(p)):
            self.R_Langevin[i] = self.k_rec*(params.N*params.N*n[i]*p[i] - self.n1*self.p1)
            
        if self.R_Langevin < 0:
            self.R_Langevin = 0
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19, 2018

@author: Timofey Golubev

This contains functions to calculate recombination rates. More types of recombination will be
added later.
"""

import numpy as np
from numba import jit

class Recombo():
    
    def __init__(self, params):
        self.R_Langevin = np.zeros(params.num_cell)
    
    @jit  
    def compute_R_Langevin(self, R_Langevin, n, p, N, k_rec, n1, p1):
        '''
        Computes bimolecular Langevin recombination rate.
        Inputs:
            R_Langevin: the empty numpy array. This is input explicitely b/c of a speedup over accessing
                        it through the recombo object.
            n: electron density
            p: hole density
            N: density of states scaling factor
            k_rec: recombination coefficient
            n1: N_LUMO*exp(-(E_LUMO - Et)/(k_B T)) number of electrons in the LUMO band when the electron’s quasi-Fermi energy
                equals the trap energy Et
            p1: N_HOMO*exp(-(Et - E_HOMO)/(k_B T)) number of holes in the HOMO band when hole’s quasi-Fermi
                energy equals Et
            n1 and p1 are defined inside of initialization.py
            
        Output: R_Langevin recombination rate array, indexed from 1.
        '''
        
        R_Langevin[1:] = k_rec*(N*N*n[1:]*p[1:] - n1*p1)
        
        # negative recombination values are unphysical
        for val in R_Langevin:
            if val < 0: 
                val = 0
                
        return R_Langevin
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19, 2018

@author: Tim
"""

import numpy as np
import constants as const
from numba import jit

class Poisson():
    
    def __init__(self, params):
        num_cell = params.num_cell
        self.epsilon = params.eps_active*np.ones(num_cell+1)  # save epsilon as attribute of the class
        self.main_diag = np.ones(num_cell)
        self.upper_diag = np.ones(num_cell-1)
        self.lower_diag = np.ones(num_cell-1) 
        
         #NOTE: PYTHON DOES RANGES AS [ )  IT DOESN'T INCLUDE THE endpoint!!!!!!!
        self.main_diag[1:num_cell] = -2*self.epsilon[1:num_cell] 
        self.upper_diag[1:num_cell-1] = self.epsilon[1:num_cell-1]  #NOTE: need to specify the indices b/c not all elements used! different array lengths!
        self.lower_diag[1:num_cell-1] = self.epsilon[1:num_cell-1] 
            
        self.rhs =  np.zeros(num_cell)
        
        self.CV = params.N*params.dx*params.dx*const.q/(const.epsilon_0*const.Vt)
               
    
    def set_rhs(self, n, p, V_left_BC, V_right_BC):
            
        self.rhs = self.CV * (n - p)
                
        self.rhs[1] -= self.epsilon[0]*V_left_BC
        self.rhs[len(self.rhs)-1] -= self.epsilon[len(self.rhs)] * V_right_BC

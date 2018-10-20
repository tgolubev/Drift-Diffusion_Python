# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19, 2018

@author: Tim
"""

import numpy as np

class Poisson():
    
    def __init(self, params):
        num_cell = params.num_cell
        self.epsilon = params.epsilon*np.ones(num_cell)  # save epsilon as attribute of the class
        self.main_diag = -2*self.epsilon*np.ones(num_cell)
        self.upper_diag = self.epsilon*np.ones(num_cell-1)
        self.lower_diag = self.epsilon*np.ones(num_cell-1)  
        self.rhs =  np.zeros(num_cell)
        
        self.CV = params.N*params.dx*params.dx*q/(epsilon_0*Vt)
        
  
    def set_rhs(self, n, p, V_left_BC, V_right_BC):
            
        for i in range(1, len(self.rhs)):
            self.rhs = self.CV * (n[i] - p[i])
                
        self.rhs[1] -= self.epsilon*V_left_BC
        self.rhs[len(self.rhs)-1] -= self.epsilon[len(self.rhs)] * V_right_BC

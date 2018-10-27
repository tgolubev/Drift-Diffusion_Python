# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19, 2018

@author: Tim
"""

import numpy as np
import constants as const
from numba import jit

class Poisson():
    '''
    This class groups all values related to the Poisson equation, making it convenient 
    to access these values through an instance of the class. Initialization of an instance of Poisson
    will also set the values of the diagonals in the Poisson matrix, since they stay constant during
    the simulation.
    '''
    
    def __init__(self, params):
        num_cell = params.num_cell
        self.epsilon = params.eps_active*np.ones(num_cell+1)  # relative dielectric constant
        self.main_diag = np.ones(num_cell)
        self.upper_diag = np.ones(num_cell-1)
        self.lower_diag = np.ones(num_cell-1) 
        
        # since the values of the Poisson matrix do not change during the simulation, we initialize
        # them only once here.
        self.main_diag[1:] = -2*self.epsilon[1:num_cell] 
        self.upper_diag[1:] = self.epsilon[1:num_cell-1]  
        self.lower_diag[1:] = self.epsilon[1:num_cell-1] 
            
        self.rhs =  np.zeros(num_cell)
        
        self.CV = params.N*params.dx*params.dx*const.q/(const.epsilon_0*const.Vt)
               
        
    @jit
    def set_rhs(self, n, p, V_left_BC, V_right_BC):
        '''
        Update the right hand side of the Poisson equation. This is done in every iteration of the 
        self consistent method due to changing charge density values and applied voltage boundary conditions.
        
        Inputs:
            n: electron density
            p: hole density
            V_left_BC: left electric potential boundary condition
            V_right_BC: right electric potential boundary condition
        '''
            
        self.rhs = self.CV * (n - p)
                
        self.rhs[1] -= self.epsilon[0] * V_left_BC
        self.rhs[-1] -= self.epsilon[-1] * V_right_BC

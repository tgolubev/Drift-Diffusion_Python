# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19, 2018

@author: Timofey Golubev

This contains everything needed to set up the continuity equation for holes 
(quasi-particle corresponding to lack of an electron), using Scharfetter-Gummel discretization.
"""

import numpy as np
import constants as const
from numba import jit

class Continuity_p():
    '''
    This class groups all values related to the hole (quasi-particle corresponding to lack of an electron)
    continuity equation, making it convenient  to access these values through an instance of the class.
    '''
    
    def __init__(self, params):

        num_cell = params.num_cell
    
        # allocate the arrays 
        self.B_p1 =  np.zeros(num_cell+1)
        self.B_p2 =  np.zeros(num_cell+1)
                
        self.main_diag = np.zeros(num_cell)
        self.upper_diag = np.zeros(num_cell-1)
        self.lower_diag = np.zeros(num_cell-1)
        self.rhs =  np.zeros(num_cell)
                
        # setup the constant arrays and variables
        self.p_mob = (params.p_mob_active/params.mobil)*np.ones(num_cell+1)
        self.Cp = params.dx*params.dx/(const.Vt*params.N*params.mobil)
        self.p_leftBC = (params.N_HOMO*np.exp(-params.phi_a/const.Vt))/params.N
        self.p_rightBC = (params.N_HOMO*np.exp(-(params.E_gap - params.phi_c)/const.Vt))/params.N
            
    @jit
    def setup_eqn(self, V, Up):
        '''
        Sets up the left and right side of the continuity matrix equation for holes. The tridiagonal matrix
        is stored in an efficient way by only storing the 3 diagonals.
        '''
        
        # update values of B_p1(V) and B_p2(V), needed in the Scharfetter-Gummel discretization
        bernoulli_fnc_p(V, self.B_p1, self.B_p2)     
        
        # set rhs
        self.rhs = -self.Cp * Up
        self.rhs[1] -= self.p_mob[0]*self.B_p1[1]*self.p_leftBC
        self.rhs[-1] -= self.p_mob[-1]*self.B_p2[-1]*self.p_rightBC
        
        # set main diagonal
        self.main_diag[1:] = -(self.p_mob[1:-1]*self.B_p2[1:-1] + self.p_mob[2:]*self.B_p1[2:])  #[1:-1] means go to 1 before last element
        
        # set lower diagonal
        self.lower_diag[1:] = self.p_mob[2:-1]*self.B_p1[2:-1]   
        
        # set upper diagonal
        self.upper_diag[1:] = self.p_mob[2:-1]*self.B_p2[2:-1] 
        

# this is defined outside of the class b/c is faster this way   
@jit(nopython = True) 
def bernoulli_fnc_p(V, B_p1, B_p2):
    '''
    This updates the values of B_p1(V) and B_p2(V) (attributes of Continuity_p class) which are 
    used in the Scharfetter-Gummel formalism of the continuity equation
    
    B_p1 = dV/(exp(dV)-1)
    B_p2 = -dV/(exp(-dV) -1) = B_p1 * exp(dV)
    
    No return value
    '''
    
    dV = np.empty(len(V))
        
    for i in range(1,len(V)):
        dV[i] = V[i] - V[i-1]
    
    B_p1[1:] = dV[1:]/(np.exp(dV[1:]) - 1.0)  #note: B_p's have length num_cell+1 since defined based on V
    B_p2[1:] = B_p1[1:]*np.exp(dV[1:])  
        


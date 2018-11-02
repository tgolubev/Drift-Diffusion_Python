# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19, 2018

@author: Timofey Golubev

This contains everything needed to set up the continuity equation for electrons, using 
Scharfetter-Gummel discretization.

"""

import numpy as np, math
import constants as const
from numba import jit

class Continuity_n():
    '''
    This class groups all values related to the electron continuity equations, making it convenient 
    to access these values through an instance of the class.
    '''
    
    def __init__(self, params):

        num_cell = params.num_cell
                
        # allocate the arrays
        self.B_n1 =  np.zeros(num_cell+1)
        self.B_n2 =  np.zeros(num_cell+1)
                
        self.main_diag = np.zeros(num_cell)
        self.upper_diag = np.zeros(num_cell-1)
        self.lower_diag = np.zeros(num_cell-1)
        self.rhs =  np.zeros(num_cell)           #right hand side
        
        # setup the constant arrays and variables
        self.n_mob = (params.n_mob_active/params.mobil)*np.ones(num_cell+1)        
        self.Cn = params.dx*params.dx/(const.Vt*params.N*params.mobil) #coeffient in front of rhs
        self.n_leftBC = (params.N_LUMO*math.exp(-(params.E_gap - params.phi_a)/const.Vt))/params.N #this is anode
        self.n_rightBC = (params.N_LUMO*math.exp(-params.phi_c/const.Vt))/params.N
        
 
    @jit    
    def setup_eqn(self, V, Un):
        '''
        Sets up the left and right side of the continuity matrix equation for electrons. The tridiagonal matrix
        is stored in an efficient way by only storing the 3 diagonals.
        '''
        
        # update values of B_n1(V) and B_n2(V), needed in the Scharfetter-Gummel discretization
        bernoulli_fnc_n(V, self.B_n1, self.B_n2) 
        
        # set rhs
        self.rhs = -self.Cn * Un                
        self.rhs[1] -= self.n_mob[0]*self.B_n2[1]*self.n_leftBC;
        self.rhs[-1] -= self.n_mob[-1]*self.B_n1[-1]*self.n_rightBC;
        
        # set main diagonal
        self.main_diag[1:] = -(self.n_mob[1:-1]*self.B_n1[1:-1] + self.n_mob[2:]*self.B_n2[2:])
        
        # set lower diagonal
        self.lower_diag[1:] = self.n_mob[2:-1]*self.B_n2[2:-1]
        
        # set upper diagonal
        self.upper_diag[1:] = self.n_mob[2:-1]*self.B_n1[2:-1]
        

@jit(nopython = True)
def bernoulli_fnc_n(V, B_n1, B_n2):
    '''
    This updates the values of B_n1(V) and B_n2(V) (attributes of Continuity_n class) which are 
    used in the Scharfetter-Gummel formalism of the continuity equation
    
    B_n1 = dV/(exp(dV)-1)
    B_n2 = -dV/(exp(-dV) -1) = B_n1 * exp(dV)
    
    No return value
    '''   
    
    dV = np.empty(len(V))
        
    for i in range(1,len(V)):
        dV[i] = V[i] - V[i-1]
        
    B_n1[1:] = dV[1:]/(np.exp(dV[1:]) - 1.0)
    B_n2[1:] = B_n1[1:]*np.exp(dV[1:])
            
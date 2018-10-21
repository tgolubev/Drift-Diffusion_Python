# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19, 2018

@author: Tim
"""

import numpy as np, math
import constants as const
from numba import jit

class Continuity_n():
    
    def __init__(self, params):

        num_cell = params.num_cell
                
        # setup the arrays
        self.n_mob = (params.n_mob_active/params.mobil)*np.ones(num_cell+1)
        self.B_n1 =  np.zeros(num_cell+1)
        self.B_n2 =  np.zeros(num_cell+1)
                
        self.main_diag = np.zeros(num_cell)
        self.upper_diag = np.zeros(num_cell-1)
        self.lower_diag = np.zeros(num_cell-1)
        self.rhs =  np.zeros(num_cell)
                
        self.Cn = params.dx*params.dx/(const.Vt*params.N*params.mobil)
        self.n_leftBC = (params.N_LUMO*math.exp(-(params.E_gap - params.phi_a)/const.Vt))/params.N #this is anode
        self.n_rightBC = (params.N_LUMO*math.exp(-params.phi_c/const.Vt))/params.N
        

    # ----------------------------------------------------------------------------   
    @jit
    def set_main_diag(self, main_diag, n_mob, B_n1, B_n2):
            
        for i in range(1, len(main_diag)):
            main_diag[i] = -(n_mob[i]*B_n1[i] + n_mob[i+1]*B_n2[i+1])
    
    @jit      
    def set_upper_diag(self, upper_diag, n_mob, B_n1):
        
        for i in range(1, len(upper_diag)):
            upper_diag[i] = n_mob[i+1]*B_n1[i+1]
    
    @jit    
    def set_lower_diag(self, lower_diag, n_mob, B_n2):

        for i in range(1, len(lower_diag)):
            lower_diag[i] = n_mob[i+1]*B_n2[i+1]
    
    @jit
    def set_rhs(self, Un):
        
        self.rhs = -self.Cn * Un
        
        #for i in range(1, len(self.rhs)):
            #self.rhs[i] = -self.Cn * Un[i]
                
        self.rhs[1] -= self.n_mob[0]*self.B_n2[1]*self.n_leftBC;
        self.rhs[len(self.rhs)-1] -= self.n_mob[len(self.rhs)]*self.B_n1[len(self.rhs)]*self.n_rightBC;
         
    def setup_eqn(self, V, Un):
        bernoulli_fnc_n(V, self.B_n1, self.B_n2)
        self.set_main_diag(self.main_diag, self.n_mob, self.B_n1, self.B_n2)
        self.set_upper_diag(self.upper_diag, self.n_mob, self.B_n1)
        self.set_lower_diag(self.lower_diag, self.n_mob, self.B_n2)
        self.set_rhs(Un)

@jit(nopython = True)
def bernoulli_fnc_n(V, B_n1, B_n2):
    dV = np.zeros(len(V))
        
    for i in range(1,len(V)):
        dV[i] = V[i] - V[i-1]
        B_n1[i] = dV[i]/(np.exp(dV[i]) - 1.0)
        B_n2[i] = B_n1[i]*np.exp(dV[i])
            
        #note: since python passes variables like by reference, if the variables are mutatable, like
        # an array, we can mutate them here and the changes will remain outside the function!            

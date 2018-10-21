# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19, 2018

@author: Tim
"""

import numpy as np
import constants as const
from numba import jit

class Continuity_p():
    
    def __init__(self, params):

        num_cell = params.num_cell
    
        # setup the arrays
        self.p_mob = (params.p_mob_active/params.mobil)*np.ones(num_cell+1)
        self.B_p1 =  np.zeros(num_cell+1)
        self.B_p2 =  np.zeros(num_cell+1)
                
        self.main_diag = np.zeros(num_cell)
        self.upper_diag = np.zeros(num_cell-1)
        self.lower_diag = np.zeros(num_cell-1)
        self.rhs =  np.zeros(num_cell)
                
        self.Cp = params.dx*params.dx/(const.Vt*params.N*params.mobil)
        self.p_leftBC = (params.N_HOMO*np.exp(-params.phi_a/const.Vt))/params.N
        self.p_rightBC = (params.N_HOMO*np.exp(-(params.E_gap - params.phi_c)/const.Vt))/params.N
            
    
    def setup_eqn(self, V, Up):
        self.bernoulli_fnc_p(V)
        self.set_main_diag()
        self.set_upper_diag()
        self.set_lower_diag()
        self.set_rhs(Up)
    
    
    # ----------------------------------------------------------------------------   
    def set_main_diag(self):
            
        for i in range(1, len(self.main_diag)):
            self.main_diag[i] = -(self.p_mob[i]*self.B_p2[i] + self.p_mob[i+1]*self.B_p1[i+1])
            
    def set_upper_diag(self):
        
        for i in range(1, len(self.upper_diag)):
            self.upper_diag[i] = self.p_mob[i+1]*self.B_p2[i+1]
            
    def set_lower_diag(self):
        
        for i in range(1, len(self.lower_diag)):
            self.lower_diag[i] = self.p_mob[i+1]*self.B_p1[i+1]
      
    def set_rhs(self, Up):
        
        for i in range(1, len(self.rhs)):
            self.rhs[i] = -self.Cp * Up[i]
                
        self.rhs[1] -= self.p_mob[0]*self.B_p1[1]*self.p_leftBC;
        self.rhs[len(self.rhs)-1] -= self.p_mob[len(self.rhs)]*self.B_p2[len(self.rhs)]*self.p_rightBC;
            
    def bernoulli_fnc_p(self, V):
        dV = np.zeros(len(V))
        
        for i in range(1,len(V)):
            dV[i] = V[i] - V[i-1]
            self.B_p1[i] = dV[i]/(np.exp(dV[i]) - 1.0)
            self.B_p2[i] = self.B_p1[i]*np.exp(dV[i])
            
            

# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19, 2018

@author: Tim
"""

import numpy as np, math
import constants as const

class Continuity_n():
    
    def __init__(self, params):

        num_cell = params.num_cell
                
        # setup the arrays
        self.n_mob = np.zeros(num_cell+1)
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
    def set_main_diag(self):
            
        for i in range(1, len(self.main_diag)):
            self.main_diag[i] = -(self.n_mob[i]*self.B_n1[i] + self.n_mob[i+1]*self.B_n2[i+1])
            
    def set_upper_diag(self):
        
        for i in range(1, len(self.upper_diag)):
            self.upper_diag[i] = self.n_mob[i+1]*self.B_n1[i+1]
            
    def set_lower_diag(self):

        for i in range(1, len(self.lower_diag)):
            self.lower_diag[i] = self.n_mob[i+1]*self.B_n2[i+1]
      
    def set_rhs(self, Un):
            
        for i in range(1, len(self.rhs)):
            self.rhs[i] = -self.Cn * Un[i]
                
        self.rhs[1] -= self.n_mob[0]*self.B_n2[1]*self.n_leftBC;
        self.rhs[len(self.rhs)-1] -= self.n_mob[len(self.rhs)]*self.B_n1[len(self.rhs)]*self.n_rightBC;
        
        
    def bernoulli_fnc_n(self, V):
        dV = np.zeros(len(V))
        
        for i in range(1,len(V)):
            dV[i] = V[i] - V[i-1]
            self.B_n1[i] = dV[i]/(np.exp(dV[i]) - 1.0)
            self.B_n2[i] = self.B_n1[i]*np.exp(dV[i])
            
            
    def setup_eqn(self, V, Un):
        self.bernoulli_fnc_n(V)
        self.set_main_diag()
        self.set_upper_diag()
        self.set_lower_diag()
        self.set_rhs(Un)


            

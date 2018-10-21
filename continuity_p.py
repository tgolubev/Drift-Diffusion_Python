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
        bernoulli_fnc_p(V, self.B_p1, self.B_p2)
        self.set_main_diag(self.main_diag,self.p_mob, self.B_p1, self.B_p2)
        self.set_upper_diag(self.upper_diag, self.p_mob, self.B_p2)
        self.set_lower_diag(self.lower_diag, self.p_mob, self.B_p1)
        self.set_rhs(Up)
    
    
    # ----------------------------------------------------------------------------   

    @jit  
    def set_rhs(self, Up):
        
        self.rhs = -self.Cp * Up
        #for i in range(1, len(self.rhs)):
            #self.rhs[i] = -self.Cp * Up[i]
                
        self.rhs[1] -= self.p_mob[0]*self.B_p1[1]*self.p_leftBC
        self.rhs[len(self.rhs)-1] -= self.p_mob[len(self.rhs)]*self.B_p2[len(self.rhs)]*self.p_rightBC
        
    @jit      
    def set_main_diag(self, main_diag, p_mob, B_p1, B_p2):
                
        tmp1 = p_mob*B_p2
        tmp2 = p_mob*B_p1
        main_diag[1:len(main_diag)] = -(tmp1[1:len(main_diag)] + tmp2[2:len(main_diag)+1])
        #for i in range(1, len(main_diag)):
         #   main_diag[i] = -(tmp1[i] + tmp2[i+1]) #-(p_mob[i]*B_p2[i] + p_mob[i+1]*B_p1[i+1])
        
    @jit   
    def set_upper_diag(self, upper_diag, p_mob, B_p2):
            
        tmp = p_mob*B_p2
        upper_diag[1:len(upper_diag)] = tmp[2:len(upper_diag)+1]  # recall python ranges omit last element
        #for i in range(1, len(upper_diag)):
        #upper_diag[i] = p_mob[i+1]*B_p2[i+1]
        
    @jit        
    def set_lower_diag(self, lower_diag, p_mob, B_p1):
            
        tmp = p_mob*B_p1
        lower_diag[1:len(lower_diag)] = tmp[2:len(lower_diag)+1]
        #for i in range(1, len(lower_diag)):
        #lower_diag[i] = p_mob[i+1]*B_p1[i+1]
    
    
# this is defined outside of the class b/c is faster this way   
@jit(nopython = True) 
def bernoulli_fnc_p(V, B_p1, B_p2):
    dV = np.zeros(len(V))
        
    for i in range(1,len(V)):
        dV[i] = V[i] - V[i-1]
        B_p1[i] = dV[i]/(np.exp(dV[i]) - 1.0)
        B_p2[i] = B_p1[i]*np.exp(dV[i])
        

 

            
            

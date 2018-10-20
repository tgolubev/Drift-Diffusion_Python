# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19, 2018

@author: Tim
"""

import numpy as np
import params

num_cell = params.num_cell
        
# setup the arrays
p_mob = np.zeros(num_cell+1)
B_p1 =  np.zeros(num_cell+1)
B_p2 =  np.zeros(num_cell+1)
        
main_diag = np.zeros(num_cell)
upper_diag = np.zeros(num_cell-1)
lower_diag = np.zeros(num_cell-1)
rhs =  np.zeros(num_cell)
        
Cp = params.dx*params.dx/(Vt*params.N*params.mobil)
p_leftBC = (params.N_HOMO*exp(-params.phi_a/Vt))/params.N
p_rightBC = (params.N_HOMO*exp(-(params.E_gap - params.phi_c)/Vt))/params.N
        

def setup_eqn(V, Up):
    bernoulli_fnc_p(V)
    set_main_diag()
    set_upper_diag()
    set_lower_diag()
    set_rhs(Up)


# ----------------------------------------------------------------------------   
def set_main_diag():
        
    for i in range(1, len(main_diag)):
        main_diag[i] = -(p_mob[i]*B_p2[i] + p_mob[i+1]*B_p1[i+1])
        
def set_upper_diag():
    
    for i in range(1, len(upper_diag)):
        upper_diag[i] = p_mob[i+1]*B_p2[i+1]
        
def set_lower_diag():
    
    for i in range(1, len(lower_diag)):
        lower_diag[i] = p_mob[i+1]*B_p1[i+1]
  
def set_rhs(Un):
        
    for i in range(1, len(rhs)):
        rhs = -Cp * Up[i]
            
    rhs[1] -= p_mob[0]*B_p1[1]*p_leftBC;
    rhs[len(rhs)-1] -= p_mob[rhs.size()]*B_p2[rhs.size()]*p_rightBC;
    
    
def bernoulli_fnc_p(V):
    dV = np.zeros(len(V))
    
    for i in range(1,len(V)):
        dV[i] = V[i] - V[i-1]
        B_p1[i] = dV[i]/(np.exp(dV[i]) - 1.0)
        B_p2[i] = B_p1[i]*np.exp(dV[i])
        
        

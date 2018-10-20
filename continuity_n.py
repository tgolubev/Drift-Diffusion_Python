# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19, 2018

@author: Tim
"""

import numpy as np
import params

num_cell = params.num_cell
        
# setup the arrays
n_mob = np.zeros(num_cell+1)
B_n1 =  np.zeros(num_cell+1)
B_n2 =  np.zeros(num_cell+1)
        
main_diag = np.zeros(num_cell)
upper_diag = np.zeros(num_cell-1)
lower_diag = np.zeros(num_cell-1)
rhs =  np.zeros(num_cell)
        
Cn = params.dx*params.dx/(Vt*params.N*params.mobil)
n_leftBC = (params.N_LUMO*exp(-(params.E_gap - params.phi_a)/Vt))/params.N #this is anode
n_rightBC = (params.N_LUMO*exp(-params.phi_c/Vt))/params.N
        

def setup_eqn(V, Un):
    bernoulli_fnc_n(V)
    set_main_diag()
    set_upper_diag()
    set_lower_diag()
    set_rhs(Un)


# ----------------------------------------------------------------------------   
def set_main_diag():
        
    for i in range(1, len(main_diag)):
        main_diag[i] = -(n_mob[i]*B_n1[i] + n_mob[i+1]*B_n2[i+1])
        
def set_upper_diag():
    
    for i in range(1, len(upper_diag)):
        upper_diag[i] = n_mob[i+1]*B_n1[i+1]
        
def set_lower_diag():
    
    for i in range(1, len(lower_diag)):
        lower_diag[i] = n_mob[i+1]*B_n2[i+1]
  
def set_rhs(Un):
        
    for i in range(1, len(rhs)):
        rhs = -Cn * Un[i]
            
    rhs[1] -= n_mob[0]*B_n2[1]*n_leftBC;
    rhs[len(rhs)-1] -= n_mob[len(rhs)]*B_n1[len(rhs)]*n_rightBC;
    
    
def bernoulli_fnc_n(V):
    dV = np.zeros(len(V))
    
    for i in range(1,len(V)):
        dV[i] = V[i] - V[i-1]
        B_n1[i] = dV[i]/(np.exp(dV[i]) - 1.0)
        B_n2[i] = B_n1[i]*np.exp(dV[i])
        
        

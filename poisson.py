# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19, 2018

@author: Tim
"""

import numpy as np
import params  #import another module

# Note: creating a class isn't really necessary, since we would have only 1 instance
# of that class, so just use a module.

num_cell = params.num_cell
epsilon = params.epsilon*np.ones(num_cell)  # save epsilon as attribute of the class
main_diag = -2*epsilon*np.ones(num_cell)
upper_diag = epsilon*np.ones(num_cell-1)
lower_diag = epsilon*np.ones(num_cell-1)  
rhs =  np.zeros(num_cell)
        
CV = params.N*params.dx*params.dx*q/(epsilon_0*Vt)
        
  
def set_rhs(n, p, V_left_BC, V_right_BC):
        
    for i in range(1, len(rhs)):
        rhs = CV * (n[i] - p[i])
            
    rhs[1] -= epsilon*V_left_BC
    rhs[len(rhs)-1] -= epsilon[len(rhs)] * V_right_BC

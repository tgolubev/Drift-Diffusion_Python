# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19, 2018

@author: Tim
"""
import numpy as np

def thomas_solve(system):
    
    num_elements = len(system.main_diag) - 1
    diagonal = system.main_diag  # make a copy so don't change the original
    
    x = np.zeros(num_elements+2)  # array for the solution
    
   # Forward substitution
    for i in range(2, num_elements + 1): # recall range uses [ ) ]
        cdiag_ratio = system.lower_diag[i-1]/diagonal[i-1]
        diagonal[i] -= cdiag_ratio * system.upper_diag[i-1]
        system.rhs[i] -= cdiag_ratio * system.rhs[i-1]
        
    # Backward substitution
    x[num_elements] = system.rhs[num_elements]/diagonal[num_elements] #lin eqn. corresponding to last row
    for i in range(num_elements, 1, -1):  # 3rd argument is the step (iterate down)
        x[i-1] = (system.rhs[i-1] - x[i]*system.upper_diag[i-1])/system.main_diag[i-1]
        
    return x
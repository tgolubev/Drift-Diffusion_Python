# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19, 2018

@author: Tim
"""
import numpy as np

def thomasSolve(matrix, rhs):
    
    num_elements = matrix.num_elements
    
    x = np.zeros(num_elements+2)  # array for the solution
    
   # Forward substitution
    for i in range(2, num_elements + 1): # recall range uses [ ) ]
        cdiag_ratio = matrix.lower_diag[i-1]/matrix.main_diag[i-1]
        matrix.main_diag[i] -= cdiag_ratio * matrix.upper_diag[i-1]
        rhs[i] -= cdiag_ratio * rhs[i-1]
        
    # Backward substitution
    x[num_elements] = rhs[num_elements]/matrix.main_diag[num_elements] #lin eqn. corresponding to last row
    for i in range(num_elements, 1, -1):  # 3rd argument is the step (iterate down)
        x[i-1] = (rhs[i-1] - x[i]*matrix.upper_diag[i-1])/matrix.main_diag[i-1]
        
    return x
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19, 2018

@author: Timofey Golubev

This contains an implementation of the Thomas algorithm for solving a tridiagonal matrix equation.
"""
import numpy as np
import time
from numba import jit

@jit(nopython=True, parallel = True)  
def thomas_solve(diagonal, upper, lower, rhs):
    '''
    Solves a tridiagonal matrix equation using the Thomas algorithm [1].
    Inputs: The matrix is passed in terms of 3 NumPy arrays corresponding to the upper, lower, and main
        diagonal. 
        rhs: array for the right hand side
        

    Reference:
        [1] https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    '''
    
    num_elements = len(diagonal) - 1

    #start = time.clock()
    
    diagonal = np.copy(diagonal)  #do a copy  to ensure the original diag is not changed
     
    x = np.empty(num_elements+2) 

    # Forward substitution  
    for i in range(2, num_elements + 1): # recall range uses [ ) ]
        cdiag_ratio = lower[i-1]/diagonal[i-1]
        diagonal[i] -= cdiag_ratio * upper[i-1]
        rhs[i] -= cdiag_ratio * rhs[i-1]
        
    # Backward substitution
    x[num_elements] = rhs[num_elements]/diagonal[num_elements] #lin eqn. corresponding to last row
    for i in range(num_elements, 1, -1): 
        x[i-1] = (rhs[i-1] - x[i]*upper[i-1])/diagonal[i-1]
        
    #end = time.clock() 
        
    return x
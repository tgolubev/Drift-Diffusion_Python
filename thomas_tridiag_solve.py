# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19, 2018

@author: Tim
"""
import numpy as np
import copy  # to be able to copy the main diagonal, so don't change the original

import time
import numba 
from numba import jit

# nopython=True means an error will be raised
# if fast compilation is not possible.

# no python mode is the fast mode!!, if it can't compile for no python mod, 
# numba will throw an exception!!
# right now it can't compile this function with no python...

@jit(nopython=True)  #this is a decorator which will make a compiled version of the function for speedup
def thomas_solve(diagonal, upper, lower, rhs):
    
    num_elements = len(diagonal) - 1
    
    # NOTE: in python, assignment statements do not copy objects!!!
    # if just asign a variable to another variable, it is like just creating another reference to it, 
    # and changing it, will change the original variable!!
    
    # use copy module to make copies.
    
    #start = time.clock()
    
    # can also do a numpy copy!!, instead of a built-in python copy
    diagonal = np.copy(diagonal)
    #upper = np.copy(upper)
    #lower = np.copy(lower)
    #rhs = np.copy(rhs)
    #seems cpu time is the same with either copy approach
     
    #diagonal = copy.deepcopy(system.main_diag)  # make a copy so don't change the original
    #upper = copy.deepcopy(system.upper_diag)
    #lower = copy.deepcopy(system.lower_diag)
    #rhs = copy.deepcopy(system.rhs)
    
    # NOTE: using deep copies here, instead of reaching into system., gave a tiny, improvement
    # in time, but not much.
     
    x = np.empty(num_elements+2)  # array for the solution. allocating "empty" instead of "zeros" is slightly faster

    # Forward substitution
      
    for i in range(2, num_elements + 1): # recall range uses [ ) ]
        cdiag_ratio = lower[i-1]/diagonal[i-1]
        diagonal[i] -= cdiag_ratio * upper[i-1]
        rhs[i] -= cdiag_ratio * rhs[i-1]
        
    # Backward substitution
    x[num_elements] = rhs[num_elements]/diagonal[num_elements] #lin eqn. corresponding to last row
    for i in range(num_elements, 1, -1):  # 3rd argument is the step (iterate down)
        x[i-1] = (rhs[i-1] - x[i]*upper[i-1])/diagonal[i-1]
        
    #end = time.clock() 
        
    return x
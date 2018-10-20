# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19, 2018

@author: Tim
"""

import continuity_n, continuity_p, initialization, photogeneration, poisson, recombination
import thomas_tridiag_solve, utilities

import numpy as np, time, matplotlib, math

params = initialization.Params()

num_cell = params.num_cell
Vbi = params.WF_anode - params.WF_cathode +params.phi_a +params.phi_c
num_V = math.floor((params.Va_max-params.Va_min)/params.increment) + 1

params.tolerance_eq = 100*params.tolerance_i  #note: we didn't have to declare a tolerance_eq to define it here!, even though is a params attribute

JV = open("JV.txt", "a")
JV.write("x                , y                          \n")

# --------------------------------------------------------------------------
# Construct objects
poisson = poisson.Poisson(params)
recombo = recombination.Recombo(params)
continuity_p = continuity_p.Continuity_p(params)
continuity_n = continuity_n.Continuity_n(params)






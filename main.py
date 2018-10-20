# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19, 2018

@author: Tim
"""

import continuity_n, continuity_p, initialization, photogeneration, poisson, recombination
import thomas_tridiag_solve, utilities

import numpy as np, time, matplotlib

params = initialization.Params()
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19, 2018

@author: Tim
"""

import continuity_n, continuity_p, initialization, photogeneration, poisson, recombination
import thomas_tridiag_solve as thomas, utilities, constants as const, time

import numpy as np, matplotlib, math

params = initialization.Params()

num_cell = params.num_cell
Vbi = params.WF_anode - params.WF_cathode +params.phi_a +params.phi_c
num_V = math.floor((params.Va_max-params.Va_min)/params.increment) + 1

params.tolerance_eq = 100*params.tolerance_i  #note: we didn't have to declare a tolerance_eq to define it here!, even though is a params attribute

JV = open("JV.txt", "a")
JV.write("x                , y                          \n")

# -------------------------------------------------------------------------------------------------
# Construct objects
poiss = poisson.Poisson(params)  # BE CAREFUL, DON'T CALL THE OBJECTS THE SAME NAMES AS THE MODULES!
recombo = recombination.Recombo(params)
cont_p = continuity_p.Continuity_p(params)
cont_n = continuity_n.Continuity_n(params)
photogen = photogeneration.get_photogeneration(params)

# initialize arrays
oldp = np.zeros(num_cell); newp = np.zeros(num_cell); oldn = np.zeros(num_cell);
newn = np.zeros(num_cell); oldV = np.zeros(num_cell+1); newV = np.zeros(num_cell+1); V = np.zeros(num_cell+1); 
Un = np.zeros(num_cell); Up = np.zeros(num_cell); R_Langevin = np.zeros(num_cell); photogen_rate = np.zeros(num_cell); 
Jp = np.zeros(num_cell); Jn = np.zeros(num_cell); J_total = np.zeros(num_cell); error_np_vector = np.zeros(num_cell); 

# Initial conditions
min_dense = min(cont_n.n_leftBC, cont_p.p_rightBC)
n = min_dense * np.ones(num_cell)
p = min_dense * np.ones(num_cell)

V_leftBC = -((Vbi)/(2*const.Vt) - params.phi_a/const.Vt)
V_rightBC = (Vbi)/(2*const.Vt) - params.phi_c/const.Vt
diff = (V_rightBC - V_leftBC)/num_cell
V[0] = V_leftBC  #fill V(0) here for use in Beroulli later
for i in range(1, num_cell):
    V[i] = V[i-1] + diff
V[num_cell] = V_rightBC
# is correct up to here

# note: poisson matrix is already set up when we constructed the poisson object

start =  time.time()

############################## MAIN LOOP ##########################################################

for Va_cnt in range(0, num_V + 2):
    not_converged = False 
    not_cnv_cnt = 0
    
    # equilibrium run
    if Va_cnt == 0:
        params.use_tolerance_eq()
        params.use_w_eq()
        Va = 0
    else:
        Va = params.Va_min + params.increment * (Va_cnt -1)
    
    if params.tolerance > 1e-5:
        print("Warning: Tolerance has been increased to > 1e-5. Results will be inaccurate")
        
    if Va_cnt == 1:
        params.use_tolerance_i();  #reset tolerance back
        params.use_w_i();
        photoge_rate = photogeneration.get_photogeneration(params);
    
    # Apply the voltage boundary conditions
    V_leftBC = -((Vbi-Va)/(2*const.Vt) - params.phi_a/const.Vt)
    V_rightBC = (Vbi-Va)/(2*const.Vt) - params.phi_c/const.Vt
    V[0] = V_leftBC
    V[num_cell] = V_rightBC
    
    
    error_np = 1.0
    iter = 0
    while error_np > params.tolerance:
        print(error_np)
        
        #--------------- Solve Poisson Equation---------------------------------------------------
        
        poiss.set_rhs(n, p, V_leftBC, V_rightBC) # update the rhs
        oldV = V
        newV = thomas.thomas_solve(poiss) 
        
        newV[0] = V[0]
        newV[num_cell] = V[num_cell]
        
        # Mix old and new solutions for V (for interior elements)
        if iter > 0:
             #NOTE: PYTHON DOES RANGES AS [)  IT DOESN'T INCLUDE THE endpoint!!!!!!!
            V[1:num_cell] = newV[1:num_cell]*params.w + oldV[1:num_cell]*(1.0 - params.w)
        else:
            V = newV
        
        # reset BC's
        V[0] = V_leftBC
        V[num_cell] = V_rightBC
        
        #-----------------Calculate net generation rate--------------------------------------------
        
        R_Langevin = recombo.compute_R_Langevin(params, n, p)
         #NOTE: PYTHON DOES RANGES AS [)  IT DOESN'T INCLUDE THE endpoint!!!!!!!
        Un[1:num_cell] = photogen_rate[1:num_cell] - R_Langevin[1:num_cell]
        Up[1:num_cell] = photogen_rate[1:num_cell] - R_Langevin[1:num_cell]
        
        #-----------------Solve equations for n and p----------------------------------------------
        
        cont_n.setup_eqn(V, Un)  
        oldn = n
        newn = thomas.thomas_solve(cont_n)
        
        cont_p.setup_eqn(V, Up)
        oldp = p
        newp = thomas.thomas_solve(cont_p)
        
        # if get negative p's or n's set them = 0
        for val in newp:
            if val < 0.0:
                val = 0
        for val in newn:
            if val < 0.0:
                val = 0
                
        # calculate the error
        old_error = error_np
        for i in range(1, num_cell):
            if (newp[i] != 0) and (newn[i] != 0):
                error_np_vector[i] = (abs(newp[i]-oldp[i]) + abs(newn[i]-oldn[i]))/abs(oldp[i]+oldn[i])
                
        error_np = max(error_np_vector)
        error_np_vector = np.zeros(num_cell)  # refill with 0's so have fresh one for next iter
        
        # auto decrease w if not converging
        if error_np >= old_error:
            not_cnv_cnt = not_cnv_cnt+1
        if not_cnv_cnt > 2000:
            params.reduce_w()
            params.relax_tolerance()
            
        #NOTE: PYTHON DOES RANGES AS [)  IT DOESN'T INCLUDE THE endpoint!!!!!!!
        p[1:num_cell] = newp[1:num_cell]*params.w + oldp[1:num_cell]*(1.0 - params.w)
        n[1:num_cell] = newn[1:num_cell]*params.w + oldn[1:num_cell]*(1.0 - params.w)
        p[0] = cont_p.p_leftBC
        n[0] = cont_n.n_leftBC
        # note: we are not including the right boundary point in p and n
        
        iter += 1
        
        # END of while loop
        
    # ------------- Calculate currents using Scharfetter-Gummel definition----------------------
        
    for i in range(1, num_cell):
        Jp[i] = (-(const.q*const.Vt*params.N*params.mobil/params.dx) * cont_p.p_mob[i] 
                *(p[i]*cont_p.B_p2[i] - p[i-1]*cont_p.B_p1[i]))
                    
        Jn[i] =  ((const.q*const.Vt*params.N*params.mobil/params.dx) * cont_n.n_mob[i] 
                 *(n[i]*cont_n.B_n1[i] - n[i-1]*cont_n.B_n2[i]))
                    
        J_total[i] = Jp[i] + Jn[i];
            
    #----------------------------Write to file-------------------------------------------------
    #if Va_cnt > 0:
     #   np.savetxt("JV.txt", Va)#, J_total[math.floor(params.num_cell/2)], iter)
            
    
endtime = time.time()

print(f"Total CPU time: {endtime-start}")

        
            
                






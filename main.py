# -*- coding: utf-8 -*-
"""
===================================================================================================
  Solving 1D Poisson + Drift Diffusion semiconductor equations for a solar cell using
                      Scharfetter-Gummel discretization

                        Created on Fri Oct 19, 2018

                          @author: Timofey Golubev

    The code as is will simulate a current-vs-voltage curve
     of a generic solar cell made of an active layer and electrodes.
    More equations for carrier recombination can be easily added.
    
    Photogeneration rate will be inputed from an input file (the name of the file can be specified in
    the parameters input file parameters.inp). For example the output of an optical model or an analytic 
    expression for photogeneration rate can be used. Generation rate file
    should contain num_cell-2 number of entries in a single column, corresponding to
    the the generation rate at each mesh point (except the endpoints).
    
    The code can also be applied to non-illuminated devices by
    setting photogeneration rate to 0.

===================================================================================================
"""

import continuity_n, continuity_p, initialization, photogeneration, poisson, recombination
import thomas_tridiag_solve as thomas, utilities, constants as const, time

import numpy as np, matplotlib.pyplot as plt, math

params = initialization.Params()
num_cell = params.num_cell
Vbi = params.WF_anode - params.WF_cathode +params.phi_a +params.phi_c
num_V = math.floor((params.Va_max-params.Va_min)/params.increment) + 1
params.tolerance_eq = 100*params.tolerance_i  

JV = open("JV.txt", "w") 
JV.write("# Voltage (V) \t  Current (A/m^2) \n")

# -------------------------------------------------------------------------------------------------
# Construct objects
poiss = poisson.Poisson(params)  
recombo = recombination.Recombo(params)
cont_p = continuity_p.Continuity_p(params)
cont_n = continuity_n.Continuity_n(params)
photogen = photogeneration.get_photogeneration(params)

# initialize arrays
oldp = np.zeros(num_cell); newp = np.zeros(num_cell); 
oldn = np.zeros(num_cell); newn = np.zeros(num_cell); 
oldV = np.zeros(num_cell+1); newV = np.zeros(num_cell+1); V = np.zeros(num_cell+1); 
Un = np.zeros(num_cell); Up = np.zeros(num_cell); 
R_Langevin = np.zeros(num_cell); photogen_rate = np.zeros(num_cell); 
Jp = np.zeros(num_cell); Jn = np.zeros(num_cell);
J_total = np.zeros(num_cell); error_np_vector = np.zeros(num_cell); 

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

# note: poisson matrix is already set up when we constructed the poisson object

start =  time.time()

###################################### MAIN LOOP ###################################################

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
        photogen_rate = photogeneration.get_photogeneration(params);
    
    # Apply the voltage boundary conditions
    V_leftBC = -((Vbi-Va)/(2*const.Vt) - params.phi_a/const.Vt)
    V_rightBC = (Vbi-Va)/(2*const.Vt) - params.phi_c/const.Vt
    V[0] = V_leftBC
    V[num_cell] = V_rightBC
    
    print(f"Va = {Va:2.2f} V")
    
    error_np = 1.0
    it = 0
    while error_np > params.tolerance:
        #print(error_np)
        
        #------------------------------ Solve Poisson Equation--------------------------------------
        
        poiss.set_rhs(n, p, V_leftBC, V_rightBC) 
        oldV = V
        newV = thomas.thomas_solve(poiss.main_diag, poiss.upper_diag, poiss.lower_diag, poiss.rhs) 
        
        newV[0] = V[0]
        newV[num_cell] = V[num_cell]
        
        # Mix old and new solutions for V (for interior elements), for stability of algorithm
        if it > 0:
            V[1:] = newV[1:]*params.w + oldV[1:]*(1.0 - params.w)
        else:
            V = newV
        
        # reset BC's
        V[0] = V_leftBC
        V[num_cell] = V_rightBC
               
        #---------------------------Calculate net generation rate-----------------------------------
        
        R_Langevin = recombo.compute_R_Langevin(R_Langevin, n, p, params.N, params.k_rec, params.n1, params.p1)
        Un[1:] = photogen_rate[1:] - R_Langevin[1:]
        Up[1:] = photogen_rate[1:] - R_Langevin[1:]
        
        #-----------------Solve equations for electron and hole density (n and p)-------------------
        
        cont_n.setup_eqn(V, Un)  
        oldn = n
        newn = thomas.thomas_solve(cont_n.main_diag, cont_n.upper_diag, cont_n.lower_diag, cont_n.rhs)
        
        cont_p.setup_eqn(V, Up)
        oldp = p
        newp = thomas.thomas_solve(cont_p.main_diag, cont_p.upper_diag, cont_p.lower_diag, cont_p.rhs)       
        
        # if get negative p's or n's set them = 0
        for val in newp:
            if val < 0.0:
                val = 0
        for val in newn:
            if val < 0.0:
                val = 0
                
        #--------------Calculate the difference (error) between prev. and current solution----------
        
        old_error = error_np
        for i in range(1, num_cell):
            if (newp[i] != 0) and (newn[i] != 0):
                error_np_vector[i] = (abs(newp[i]-oldp[i]) + abs(newn[i]-oldn[i]))/abs(oldp[i]+oldn[i])
        
        error_np = max(error_np_vector)
        error_np_vector = np.zeros(num_cell)  # refill with 0's so have fresh one for next it
         
        # auto decrease w if not converging
        if error_np >= old_error:
            not_cnv_cnt = not_cnv_cnt+1
        if not_cnv_cnt > 2000:
            params.reduce_w()
            params.relax_tolerance()
            
        # linear mixing of old and new solutions for stability
        p[1:num_cell] = newp[1:num_cell]*params.w + oldp[1:num_cell]*(1.0 - params.w)
        n[1:num_cell] = newn[1:num_cell]*params.w + oldn[1:num_cell]*(1.0 - params.w)
        p[0] = cont_p.p_leftBC
        n[0] = cont_n.n_leftBC
        # note: we are not including the right boundary point in p and n here
        
        it += 1
        
        # END of while loop
          
    # ------------- Calculate currents using Scharfetter-Gummel definition-------------------------
        
    for i in range(1, num_cell):
        Jp[i] = (-(const.q*const.Vt*params.N*params.mobil/params.dx) * cont_p.p_mob[i] 
                *(p[i]*cont_p.B_p2[i] - p[i-1]*cont_p.B_p1[i]))
                    
        Jn[i] =  ((const.q*const.Vt*params.N*params.mobil/params.dx) * cont_n.n_mob[i] 
                 *(n[i]*cont_n.B_n1[i] - n[i-1]*cont_n.B_n2[i]))
                    
        J_total[i] = Jp[i] + Jn[i];

    #----------------------------Write results to file----------------------------------------------
    if Va_cnt > 0:
        JV.write(f"{Va:2.2f} \t\t\t {J_total[math.floor(params.num_cell/2)]:4.8f} \n")
        
    # End of main loop
          

endtime = time.time()
print(f"Total CPU time: {endtime-start}") 
JV.close()

# Plot Results
V, J = np.loadtxt("JV.txt", usecols=(0,1), unpack = True)  # usecols specifies columns to use, unpack specifies to use tuple unpacking
plt.xlim(params.Va_min, params.Va_max)
plt.ylim(-250, 100)
plt.plot(V, J)
plt.xlabel('Voltage ($V$)')
plt.ylabel('Current ($A/m^2$)') # TeX markup
plt.grid(True)
plt.savefig("JV.jpg", dpi = 1200) #specify the dpi for a high resolution image


import numpy as np
import sys
sys.path.append('/Users/williamjo/Documents/LPI/Codes/plume_regolith_solver_v1/src')

import var_timestepping as timestep
import var_range_of_interest as bounds
import var_soil as soil
import var_impinged_gas as imp
import sys_output as out

from fun_impinged_gas import *
from fun_soil import *
from fun_nozzle import *
from fun_integration import *
from write_outputs import *

#nozzle.gamma = 1.4
#nozzle.Ma = 2
#
#
#A_Astar = A_over_Astar(nozzle.Ma, nozzle.gamma)
#T0_T = T0_over_T(nozzle.Ma, nozzle.gamma)
#p0_p = p0_over_p(nozzle.Ma, nozzle.gamma)
#rho0_rho = rho0_over_rho(nozzle.Ma, nozzle.gamma)
#pNS_p0 = pNS_over_p0(nozzle.Ma, nozzle.gamma)
#
#print(A_Astar)
#print(1/T0_T)
#print(1/p0_p)
#print(1/rho0_rho)
#print(pNS_p0)

compute_nozzle_exhaust()

print(nozzle.R_gas)
h_nozzle =  31.5 #m
min_altitude = 1 # height of lander

bounds_ej_ring = np.zeros([bounds.n_points_centerline, 4])

for i in range(1, bounds.n_points_centerline-1):
        
    bounds_ej_ring[i-1,0] = i-1
    bounds_ej_ring[i-1,1] = bounds.r_centerlines_array[i-1] 
    bounds_ej_ring[i-1,2] = soil.r_midpoint[i-1]
    bounds_ej_ring[i-1,3] = bounds.r_centerlines_array[i]

write_bounds(bounds_ej_ring)
    
soil.avg_rho_p_excavated = np.ones(bounds.n_points_centerline-1) * compute_soil_density_depth(0)

print(soil.avg_rho_p_excavated)
# ------------------------------------------------------------------------------------------------
# Step 3: Given a height h (input), array containing the radial distances from plume centerline (implict)
# And excavation depths (Initial Condition: h_excavated = 0 for all distances away from centerline)
# 
# Solve For:
# The mass erosion rate flux (m_dot, kg/m^2s), threshold energy (E_th), cohesion energy (Alpha) (implict)
# The mass erosion rate (kg/s) at each timestep for each differential ring (implicit)
# The instantaneous total mass eroded in each ring for one timestep, given a timestep size (dt) (implicit)
# The cumulative total mass eroded at each differential distance (implicit)
# ------------------------------------------------------------------------------------------------
compute_mass_eroded_quantities(h_nozzle)

#print(nozzle.m_dot_e) 
#print(nozzle.P_e)
#print(nozzle.T_e)
#print(nozzle.rho_e) 
#print(nozzle.v_e)
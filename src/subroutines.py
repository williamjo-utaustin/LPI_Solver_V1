import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import quad
from scipy.optimize import fsolve

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

import csv
import time

def loop():
    
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


    for t in range(0,6100):

        # ----------------------------------------------------------------------------------------------
        # Step 1: Print the current timestep 
        # ----------------------------------------------------------------------------------------------        
        print_timestep(t, h_nozzle)
        # ----------------------------------------------------------------------------------------------        

        # ----------------------------------------------------------------------------------------------
        # Step 2: Output the ejecta properties for the current timestep
        # ----------------------------------------------------------------------------------------------        
        plot_erosion_profile(t, soil.r_midpoint, soil.h_excavated_mid)
        # ----------------------------------------------------------------------------------------------        

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
        # ----------------------------------------------------------------------------------------------        
        
        # ----------------------------------------------------------------------------------------------
        # Step 4: Output the pressure, temperature, and velocity of the gas variables
        # ----------------------------------------------------------------------------------------------
        output_gas_props() 
        # ----------------------------------------------------------------------------------------------        

        # ----------------------------------------------------------------------------------------------
        # Step 5: The following block solves for the excavated height given a differential ring dA, its area,
        # and the instantaneous mass excavated from the prior timestep.
        # Note: need to turn block below into a function but it doesn't seem to be working.
        # 
        # Solve for: The excavated heights h at the midpoint location of each bin
        # ----------------------------------------------------------------------------------------------
        solve_excavation_depths()
        # ----------------------------------------------------------------------------------------------        
        
        # ----------------------------------------------------------------------------------------------
        # Step 6: The following block updates the excavation depths of each bin with new values
        # ----------------------------------------------------------------------------------------------
        set_new_excavation_depths()
        # ----------------------------------------------------------------------------------------------
        
        # ----------------------------------------------------------------------------------------------
        # Step 7: Determine the ejecta properties (Ejecta Velocity per Particle Size per Distance from Centerline) 
        # ----------------------------------------------------------------------------------------------        
        u_ej_timestep, offset_ej_dist_timestep, offset_ej_time_timestep, ej_timestep_props = determine_ejecta_props(t)
        # ----------------------------------------------------------------------------------------------        

        # ----------------------------------------------------------------------------------------------
        # Step 8: Output the ejecta properties for the next timestep
        # ----------------------------------------------------------------------------------------------        
        write_ejecta_props(t, h_nozzle, ej_timestep_props, u_ej_timestep, offset_ej_dist_timestep, offset_ej_time_timestep)
        # ----------------------------------------------------------------------------------------------        

        # ----------------------------------------------------------------------------------------------
        # Step 9: Update the new nozzle velocity
        # ----------------------------------------------------------------------------------------------
        h_nozzle = update_nozzle_height(h_nozzle)
        # ----------------------------------------------------------------------------------------------        

        if(h_nozzle < min_altitude):
            break
    return None
    
# ----------------------------------
# need to modify this later
# ----------------------------------
#def plot_initial_disturbed_altitude():
#    
#    for j in range(0,bounds.n_points_altitude):
#
#        h_nozzle = bounds.h_nozzle_array[j]
#        m_dot = np.zeros(int(bounds.n_points_centerline))
#
#        for i in range(0, bounds.n_points_centerline):
#
#            r_centerline = bounds.r_centerlines_array[i]
#            compute_impinged_gas(h_nozzle, r_centerline)
#            m_dot[i] = compute_mass_erosion_rate()
#
#        plt.plot(bounds.r_centerlines_array, m_dot, linewidth = 3, label = str(h_nozzle) + 'm')
#
#    plt.xlim(0,10000)
#    plt.ylim(0,7)
#    plt.xticks(fontsize = 18)
#    plt.yticks(fontsize = 18)
#    plt.xlabel("Distance from Centerline (m)", fontsize = 18)
#    plt.ylabel("Unit Area nozzle.Mass Flow Rate ($kg/{m^2s}$)", fontsize = 18)
#    plt.legend(title = "Start Altitude (m)", title_fontsize= 16, fontsize = 16 , loc = 'upper right')
#    plt.show()
    
# ----------------------------------------------------------------------------------
# Next Steps
# ----------------------------------------------------------------------------------
# Move the spacecraft to the next position (update the h and r values)
# Update the new energy threshold values for the next timestep in each radial bin
# Update the density of the regolith at the new depth in each radial bin
# Solve for the mass flow rate per unit area using updated values and repeat process
# Start to derive launch velocities and theoretical launch angles for particles
# Think about how you can get mass launched to number of particles launched at time t
# Given mass ejected, velocity, angle, density, particle size distribution, create Monte Carlo Code
# Needs a f_num value for each particle
# Solve for the flux per unit area, Mie Scattering, grain heating etc. etc.
# ----------------------------------------------------------------------------------

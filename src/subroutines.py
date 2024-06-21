import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import quad
from scipy.optimize import fsolve

import var_timestepping as timestep
import var_range_of_interest as bounds
import var_soil as soil
import var_impinged_gas as imp

from fun_impinged_gas import *
from fun_soil import *
from fun_nozzle import *
from fun_integration import *

import csv
import time

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



def loop():
    #d_particle = np.asarray([multiplier * magnitude for magnitude in [1E-6, 1E-5, 1E-4, 1E-3] for multiplier in [1, 2, 3, 4, 5, 6, 7, 8, 9]])

    # initialize_variables()
    # start from t = 0 or from a restart file

    # We need to check if the heights from previous timesteps are being saved onto the computer.
    # From analysis we see that this does not happen.

    # set the timestep (seconds)
    h_nozzle = 25 #m

    bounds_ej_ring = np.zeros([bounds.n_points_centerline, 4])

    for i in range(1, bounds.n_points_centerline-1):
        
        bounds_ej_ring[i-1,0] = i-1
        bounds_ej_ring[i-1,1] = bounds.r_centerlines_array[i-1] 
        bounds_ej_ring[i-1,2] = soil.r_midpoint[i-1]
        bounds_ej_ring[i-1,3] = bounds.r_centerlines_array[i]

    np.savetxt("output/ejecta_ring_bounds.csv", bounds_ej_ring, delimiter=',', fmt=' '.join(['%i'] + ['%.4e']*3), header = '# Ejecta Ring Radial Distances Away from Plume Centerline (Index, Min (m), Midpoint (m), Max (m))')
    print("Timestep,", "Depth Excavated,", "Threshold Energy,", "E_down,", "Alpha,","Mdot_flux", "M_area_eroded_inst", "Mdot_cumulative" )
    
    soil.avg_rho_p_excavated = np.ones(bounds.n_points_centerline-1) * compute_soil_density_depth(0)
    
    for t in range(0,6100):

        print("Writing at time t =", t * timestep.delta_t, "Nozzle Height: ", h_nozzle)

        # ------------------------------------------------------------------------------------------------
        #  Given a height h (input), array containing the radial distances from plume centerline (implict)
        #  And excavation depths (Initial Condition: h_excavated = 0 for all distances away from centerline)
        # 
        #  Solve For:
        #  The mass erosion rate flux (m_dot, kg/m^2s), threshold energy (E_th), cohesion energy (Alpha) (implict)
        #  The mass erosion rate (kg/s) at each timestep for each differential ring (implicit)
        #  The instantaneous total mass eroded in each ring for one timestep, given a timestep size (dt) (implicit)
        #  The cumulative total mass eroded at each differential distance (implicit)
        # ------------------------------------------------------------------------------------------------
        compute_mass_eroded_quantities(h_nozzle)
        

        # ----------------------------------------------------------------------------------------------
        # ----------------------------------------------------------------------------------------------
                
        #for i in range(0,bounds.n_points_centerline-1):
        #    print(i, soil.r_midpoint[i], imp.v_gas_arr[i], imp.p_gas_arr[i], imp.rho_gas_arr[i], imp.T_gas_arr[i])

        # ----------------------------------------------------------------------------------------------
        # The following block solves for the excavated height given a differential ring dA, its area,
        # and the instantaneous mass excavated from the prior timestep.
        # Note: need to turn block below into a function but it doesn't seem to be working.
        # 
        # Solve for: The excavated heights h at the midpoint location of each bin
        # Additional Needs: Don't integrate stuff that has not moved!
        # ----------------------------------------------------------------------------------------------
        
        solve_excavation_depths()
        
        # ----------------------------------------------------------------------------------------------
        # ----------------------------------------------------------------------------------------------

        # compute the total mass eroded in each ring for one timestep

        if(np.mod(t,1)==0):
            print("Plotting Figure at Timestep: ", str(t))
            fig = plt.figure(figsize=(10, 6))
            plt.title("Time = "+str(t * timestep.delta_t) + " s", fontsize = 18)
            plt.scatter(soil.r_midpoint, soil.h_excavated_mid, linewidth = 2)
            #plt.semilogy(soil.r_midpoint, soil.h_excavated_mid, linewidth = 2)
            plt.xlabel("Distance from Plume Centerline (m)", fontsize = 18)
            plt.ylabel("Excavation Depth (m)", fontsize = 18)
            plt.xticks(fontsize = 18)
            plt.yticks(fontsize = 18)
            #plt.xlim(0,10)
            plt.xlim(0,100)
            plt.ylim(0.25, 0)
            #plt.gca().set_aspect('equal')
            plt.savefig("output/growth_"+str(t)+".png", dpi = 300)
            plt.close()

        # ----------------------------------------------------------------------------------------------
        # The following block updates the excavation depths of each bin with new values
        # ----------------------------------------------------------------------------------------------
        set_new_excavation_depths()
        # ----------------------------------------------------------------------------------------------
        
        # ----------------------------------------------------------------------------------------------
        # Determine the ejecta properties (Ejecta Velocity per Particle Size per Distance from Centerline) 
        # ----------------------------------------------------------------------------------------------        
        u_ej_timestep, offset_ej_dist_timestep, offset_ej_time_timestep, ej_timestep_props = determine_ejecta_props(t)

        # saving values
        np.savetxt("output/ejecta_properties_"+str(t)+".csv", ej_timestep_props, delimiter=',', fmt=' '.join(['%i'] + ['%.8e']*4), header = '# header = Nozzle Height (m) ' + str(h_nozzle)+', Offset angle (3 deg), ' + 'Bin Index, Bin Midpoint Location (m), Excavated Density (kg/m^3), Instantaneous Mass Excavated (kg), Total Height Excavated (m)')
        np.savetxt("output/ejecta_velocities_"+str(t)+".csv", u_ej_timestep, delimiter=',', fmt=' '.join(['%i'] + ['%.4e']*np.size(soil.d_particle)), header = '# Nozzle Height (m) ' + str(h_nozzle)+', Offset angle (3 deg), '+' Ejecta Speeds by Bin Index and Particle Size (O(1), O(10), O(100), O(1000) microns) in sets of 1, 2, 3, 4, 5, 6, 7, 8, 9')
        np.savetxt("output/ejecta_offset_distances_"+str(t)+".csv", offset_ej_dist_timestep, delimiter=',', fmt=' '.join(['%i'] + ['%.4e']*np.size(soil.d_particle)), header = '# header = Nozzle Height (m) ' + str(h_nozzle)+', Offset angle (3 deg), ' + 'Ejecta Distance Offset by Bin Index and Particle Size (O(1), O(10), O(100), O(1000) microns) in sets of 1, 2, 3, 4, 5, 6, 7, 8, 9')
        np.savetxt("output/ejecta_offset_time_"+str(t)+".csv", offset_ej_time_timestep, delimiter=',', fmt=' '.join(['%i'] + ['%.4e']*np.size(soil.d_particle)), header = '# header = Nozzle Height (m) ' + str(h_nozzle)+', Offset angle (3 deg), ' + 'Ejecta Time Offset by Bin Index and Particle Size (O(1), O(10), O(100), O(1000) microns) in sets of 1, 2, 3, 4, 5, 6, 7, 8, 9')

        # update the new nozzle velocity
        h_nozzle = update_nozzle_height(h_nozzle)

    return None
    
    
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



    #print("yes")
    #plt.semilogy(r_midpoint, h_excavated_mid)
    #plt.xlim(bounds.min_centerline, bounds.max_centerline)
    #plt.ylim(1E-8, 1E1)
    #plt.show()
import numpy as np
from scipy.integrate import quad
from scipy.optimize import fsolve

import matplotlib.pyplot as plt

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
    d_particle = np.asarray([multiplier * magnitude for magnitude in [1E-6, 1E-5, 1E-4, 1E-3] for multiplier in [1, 2, 3, 4, 5, 6, 7, 8, 9]])

    # initialize_variables()
    # start from t = 0 or from a restart file

    # unload the restart file

    # set the timestep (seconds)
    h_nozzle = 20 #m

    bounds_ej_ring = np.zeros([bounds.n_points_centerline, 4])

    for i in range(1, bounds.n_points_centerline-1):
        bounds_ej_ring[i-1,0] = i-1
        bounds_ej_ring[i-1,1] = bounds.r_centerlines_array[i-1] 
        bounds_ej_ring[i-1,2] = soil.r_midpoint[i-1]
        bounds_ej_ring[i-1,3] = bounds.r_centerlines_array[i]

    np.savetxt("output/ejecta_ring_bounds.csv", bounds_ej_ring, delimiter=',', fmt=' '.join(['%i'] + ['%.4e']*3), header = '# Ejecta Ring Radial Distances Away from Plume Centerline (Index, Min (m), Midpoint (m), Max (m))')


    print("Timestep,", "Depth Excavated,", "Threshold Energy,", "E_down,", "Alpha,","Mdot_flux", "M_area_eroded_inst", "Mdot_cumulative" )
    
    for t in range(0,1000):

        print("Writing at time t =", t * timestep.delta_t)

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
        h_excavated_mid_old = np.zeros(bounds.n_points_centerline-1)

        for i in range(0,bounds.n_points_centerline-1):

            # save excavation depth for current timestep
            h_excavated_mid_old[i] = soil.h_excavated_mid[i]

            # solve for the surficial density (kg/m^2)
            soil.surf_den_excavated[i] = (soil.m_excavated_inst[i])/soil.ring_area[i]

            # solve for the depth z of excavation by integration and solving for the upper bound
            def func(z):
                y, err = quad(compute_soil_density_depth, soil.h_excavated_mid[i], z)
                return soil.surf_den_excavated[i] - y
            
            sol = fsolve(func, 1E-8)

            soil.h_excavated_mid[i] = sol[0]
        # ----------------------------------------------------------------------------------------------
        # ----------------------------------------------------------------------------------------------

        fig = plt.figure(figsize=(10, 10))
        plt.title("Time = "+str(t * timestep.delta_t) + " s", fontsize = 18)
        plt.plot(soil.r_midpoint, soil.h_excavated_mid, linewidth = 2)
        #plt.semilogy(soil.r_midpoint, soil.h_excavated_mid, linewidth = 2)
        plt.xlabel("Distance from Plume Centerline (m)", fontsize = 18)
        plt.ylabel("Excavation Depth (m)", fontsize = 18)
        plt.xticks(fontsize = 18)
        plt.yticks(fontsize = 18)
        plt.xlim(0,10)
        #plt.xlim(0,6000)
        plt.ylim(1E0, 1E-5)
        plt.savefig("output/growth_"+str(t)+".png", dpi = 300)
        plt.close()

        # ----------------------------------------------------------------------------------------------
        # The following block updates the excavation depths of each bin with new values
        # ----------------------------------------------------------------------------------------------
        set_new_excavation_depths()
        # ----------------------------------------------------------------------------------------------
        # ----------------------------------------------------------------------------------------------


        # ----------------------------------------------------------------------------------------------        
        # The following block solves for the Gthe average density of lunar material
        # excavated at each timestep for material affected by the plume
        # 
        # Solve for: Average excavated density of lunar material for bins excavated by
        # the plume (Should increase as dig further down, will be used for MC simulation)
        # ----------------------------------------------------------------------------------------------        
        #avg_rho_p_excavated = np.empty(bounds.n_points_centerline-1)

        if(t == 0):
            avg_rho_p_excavated = np.ones(bounds.n_points_centerline-1) * compute_soil_density_depth(0)

        avg_rho_p_excavated_old = avg_rho_p_excavated.copy()
        #print(avg_rho_p_excavated_old)

        index_count = np.zeros(bounds.n_points_centerline-1)
        diff = np.zeros(bounds.n_points_centerline-1)
        total_affected_indicies = 0

        # compute the average density of the material that has already been excavated
        for i in range(0,bounds.n_points_centerline-1):
            
            
            if((soil.h_excavated_mid[i] - h_excavated_mid_old[i]) > 0):
                avg_rho_p_excavated[i],err = quad(compute_soil_density_depth, h_excavated_mid_old[i], soil.h_excavated_mid[i]) * 1/(soil.h_excavated_mid[i] - h_excavated_mid_old[i])
            else:
                avg_rho_p_excavated[i] = compute_soil_density_depth(soil.h_excavated_mid[i])

            diff[i] = np.abs(avg_rho_p_excavated[i] - avg_rho_p_excavated_old[i])        

            if(np.abs(avg_rho_p_excavated[i] - avg_rho_p_excavated_old[i] > 0)):
                index_count[total_affected_indicies] = i
                total_affected_indicies = total_affected_indicies + 1

        #print(index_count)
        #print(total_affected_indicies)
        #print(diff)
        # ----------------------------------------------------------------------------------------------        
        # ----------------------------------------------------------------------------------------------        
        # Compute the ejecta speeds of each grain launched
        # Compute the properties of the gas at the midpoint of each bin
        rho_gas_arr_mid = np.zeros(bounds.n_points_centerline-1)
        u_gas_arr_mid = np.zeros(bounds.n_points_centerline-1)
        T_gas_arr_mid = np.zeros(bounds.n_points_centerline-1)
        P_gas_arr_mid = np.zeros(bounds.n_points_centerline-1)
    
        u_ej_timestep = np.zeros([total_affected_indicies-1, np.size(d_particle) + 1])
        offset_ej_dist_timestep = np.zeros([total_affected_indicies-1, np.size(d_particle) + 1])
        offset_ej_time_timestep = np.zeros([total_affected_indicies-1, np.size(d_particle) + 1])
        ej_timestep_props = np.zeros([total_affected_indicies-1, 5])


        for i in range(0,total_affected_indicies-1):

            i_n = int(index_count[i])
            ip1 = int(index_count[i+1])

            ej_timestep_props[i, 0] = i_n
            ej_timestep_props[i, 1] = soil.r_midpoint[i_n]
            ej_timestep_props[i, 2] = avg_rho_p_excavated[i_n]
            ej_timestep_props[i, 3] = soil.m_excavated_inst[i_n]
            ej_timestep_props[i ,4] = soil.h_excavated_mid[i_n]


            print("Integrating from index ", i_n, "to ", ip1, "Iteration #", i+1,"out of ",total_affected_indicies-1)

            rho_gas_arr_mid[i_n] = (imp.rho_gas_arr[i_n] + imp.rho_gas_arr[ip1])/2
            u_gas_arr_mid[i_n] = (imp.v_gas_arr[i_n] + imp.v_gas_arr[ip1])/2
            T_gas_arr_mid[i_n] = (imp.T_gas_arr[i_n] + imp.T_gas_arr[ip1])/2
            P_gas_arr_mid[i_n] = (imp.p_gas_arr[i_n] + imp.p_gas_arr[ip1])/2
            
            u_ej_timestep[i,0] = i_n
            offset_ej_dist_timestep[i,0] = i_n
            offset_ej_time_timestep[i,0] = i_n

            for j in range(0,np.size(d_particle)):
                
                d_p = d_particle[j]
                x_p = 0
                y_p = 0
                u_p = 0
                t_sub = 0

                # number of array points
                du = 0.01
                n_timesteps = (int(u_gas_arr_mid[i_n]/du) - 1)

                t_array = np.zeros(n_timesteps)
                u_p_array = np.zeros(n_timesteps)
    
                x_p_array = np.zeros(n_timesteps)
                y_p_array = np.zeros(n_timesteps)

                for ts in range (0, n_timesteps):

                    u_p_array[ts] = u_p
                    x_p_array[ts] = x_p
                    y_p_array[ts] = y_p    
    
                    t_array[ts] = t_sub
        
                    if (y_p > 0.03):
                        u_ej_timestep[i,j+1] = u_p
                        offset_ej_dist_timestep[i,j+1] = x_p
                        offset_ej_time_timestep[i,j+1] = t
                        break

                    delta_t = du/(dudt(u_p, ts, d_p, avg_rho_p_excavated[i_n], u_gas_arr_mid[i_n], P_gas_arr_mid[i_n], rho_gas_arr_mid[i_n], T_gas_arr_mid[i_n]))

                    u_p = u_p + du
                    t_sub = t_sub + delta_t

                    x_p = x_p + (u_p * delta_t) * np.cos(3 * np.pi/180)
                    y_p = y_p + (u_p * delta_t) * np.sin(3 * np.pi/180)
    
                
                #for k in range(0,timestep.n_sub_timesteps):
                #    t_sub = t_sub + timestep.dt 
                #    u_p_save = u_p
                #    u_p = rk_4(dudt, u_p, t, timestep.dt, d_p, P_gas_arr_mid[i_n], rho_gas_arr_mid[i_n],  T_gas_arr_mid[i_n], u_gas_arr_mid[i_n], avg_rho_p_excavated[i_n])

                #    x_n = x_n + u_p * timestep.dt 
                #    if(x_n * np.tan(3 * np.pi/180) > timestep.baseline_height or (u_p-u_p_save)/u_p < 0.00001):
                #        
                #        u_ej_timestep[i,j+1] = u_p
                #        offset_ej_dist_timestep[i,j+1] = x_n
                #        offset_ej_time_timestep[i,j+1] = t_sub
                #        break
            
            #print(i, avg_rho_p_excavated[i], rho_gas_arr_mid[i], v_gas_arr_mid[i], T_gas_arr_mid[i], v_T_arr_mid[i], P_gas_arr_mid[i])
            # only investigate the particles that have been moved

        #print(P_gas_arr_mid[10], u_gas_arr_mid[10], T_gas_arr_mid[10], rho_gas_arr_mid[10])
        
        # write output to data file

        print(t)

        # saving values
        np.savetxt("output/ejecta_properties_"+str(t)+".csv", ej_timestep_props, delimiter=',', fmt=' '.join(['%i'] + ['%.8e']*4), header = '# header = Nozzle Height (m) ' + str(h_nozzle)+', Offset angle (3 deg), ' + 'Bin Index, Bin Midpoint Location (m), Excavated Density (kg/m^3), Instantaneous Mass Excavated (kg), Total Height Excavated (m)')
        np.savetxt("output/ejecta_velocities_"+str(t)+".csv", u_ej_timestep, delimiter=',', fmt=' '.join(['%i'] + ['%.4e']*np.size(d_particle)), header = '# Nozzle Height (m) ' + str(h_nozzle)+', Offset angle (3 deg), '+' Ejecta Speeds by Bin Index and Particle Size (O(1), O(10), O(100), O(1000) microns) in sets of 1, 2, 3, 4, 5, 6, 7, 8, 9')
        np.savetxt("output/ejecta_offset_distances_"+str(t)+".csv", offset_ej_dist_timestep, delimiter=',', fmt=' '.join(['%i'] + ['%.4e']*np.size(d_particle)), header = '# header = Nozzle Height (m) ' + str(h_nozzle)+', Offset angle (3 deg), ' + 'Ejecta Distance Offset by Bin Index and Particle Size (O(1), O(10), O(100), O(1000) microns) in sets of 1, 2, 3, 4, 5, 6, 7, 8, 9')
        np.savetxt("output/ejecta_offset_time_"+str(t)+".csv", offset_ej_time_timestep, delimiter=',', fmt=' '.join(['%i'] + ['%.4e']*np.size(d_particle)), header = '# header = Nozzle Height (m) ' + str(h_nozzle)+', Offset angle (3 deg), ' + 'Ejecta Time Offset by Bin Index and Particle Size (O(1), O(10), O(100), O(1000) microns) in sets of 1, 2, 3, 4, 5, 6, 7, 8, 9')

        #with open("output/ejecta_velocities_"+str(t)+".csv","w+") as my_csv:
        #    csvWriter = csv.writer(my_csv,delimiter=',')
        #    csvWriter.writerows(u_ej_timestep)






        
        # update the new nozzle velocity
        h_nozzle = update_nozzle_height(h_nozzle)
        
        # print statisics
        #if(np.mod(t,10)==0):
        #    print(t, soil.h_excavated_bounds[860], soil.E_th[860], imp.E_downward, soil.alpha[860], soil.m_dot_area_eroded[860], soil.m_excavated_inst[860], soil.m_excavated_cumulative[860])





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
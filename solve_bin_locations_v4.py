import numpy as np
import matplotlib.pyplot as plt
import time
import sys

sys.path.append('/Users/williamjo/Documents/LPI/Codes/plume_regolith_solver_v1/src')

import var_nozzle as nozzle
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
import var_solve_ejection as sol_ej

#import var_soil as soil
#import var_constants as cs
#from fun_regolith_distributions import *
#from fun_farfieldbounds import *
#from write_outputs import *


print_dataset_1 = True
print_dataset_2 = True

#file_folder = "starship_H30_lambda1p5_m1e05"
file_folder = "starship_H60_lambda10_m1e05"
identifier = "_scaling_9"
max_timestep = 286

#file_folder = "bluemoon_H60_lambda2_m1p5e04"
#identifier = "_scaling_9"
#max_timestep = 860
#file_folder = "starship_raptor_nominal_50"
#identifier = "_scaled_dec2024_scaling_9"
#max_timestep = 867

#file_folder = "apollo_vel_slow_alt_mid"
#identifier = "_4_scaling_9"
#max_timestep = 315

#file_folder = "starship_raptor_nominal_50"
#identifier = "_scaled_scaling_9"
#max_timestep = 886


# Assume a 3 degree ejection angle (most probable) and orbital mechanics
range_data = np.genfromtxt("data/ejecta_ranges.csv", delimiter=",")
#sol_ej.ejecta_ring_bounds = np.genfromtxt("output_apollo/output_apollo_raw_data/apollo_ejecta_ring_bounds.csv", delimiter = ' ')

print(file_folder+identifier+"/"+file_folder+"_ejecta_ring_bounds.csv")
sol_ej.ejecta_ring_bounds = np.genfromtxt(file_folder+identifier+"/"+file_folder+"_ejecta_ring_bounds.csv", delimiter = ' ')
print(sol_ej.range_bounds)

#input("Hold")
if(sol_ej.plot_only is False):


    # -----------------------------------------------
    # save the range bounds files (for further post-processing)
    # -----------------------------------------------
    np.savetxt(file_folder+identifier+"/"+file_folder+"_post_processing/range_bounds.csv", sol_ej.range_bounds_mid, delimiter=",")
    # -----------------------------------------------
    # -----------------------------------------------
    

    # -----------------------------------------------            
    # loop through particle launch time
    # -----------------------------------------------            
    for t in range(0,max_timestep):

        if(np.mod(t,1)==0): 
            print("Time (s) = ", t * 0.01)


        # -------------------------------------------------------
        # Open the ejecta properties and ejecta velocities files
        # Ejecta velocities will have the size 37, with the zeroth column being a dummy index
        # Columns 1 - 36 will represent the speed of the ejection particle by size
        # Rows 0 - 19 represent the index of the launch location from the centerline
        # -------------------------------------------------------
        #ejecta_properties = np.genfromtxt("output_apollo/output_apollo_raw_data/apollo_ejecta_properties_"+str(t)+".csv", delimiter = " ")
        #ejecta_velocities = np.genfromtxt("output_apollo/output_apollo_raw_data/apollo_ejecta_velocities_"+str(t)+".csv", delimiter = " ")
        
        ejecta_properties = np.genfromtxt(file_folder+identifier+"/"+file_folder+"_ejecta_properties_"+str(t)+".csv", delimiter = " ")
        ejecta_velocities = np.genfromtxt(file_folder+identifier+"/"+file_folder+"_ejecta_velocities_"+str(t)+".csv", delimiter = " ")
        # -------------------------------------------------------
        
        #print(ejecta_properties)
        #print(ejecta_velocities)

        #input("Pause")

        # ---------------------------------------------------------------------------------
        # extract the instantaneous mass excavated data for each radial ring at time i * dt
        # In ejecta_properties, the file format is as follows, where (" ") represents the index of the column
        # Bin Index(0), Bin Midpoint Location (1), Excavated Density (2), Instantaneous Mass Excavated (3), Total Height Excavated (4)
        # ---------------------------------------------------------------------------------
        ejecta_indices = ejecta_properties[:,0]
        ejecta_mass_excavated = ejecta_properties[:,3]
        ejecta_rho_excavated = ejecta_properties[:,2] 

        #print(ejecta_indices)
        #print(ejecta_mass_excavated)
        #print(ejecta_rho_excavated)

        # convert mass to volume by dividing by the density
        # get volume of each particle size (1-9), (10 - 90), ...
        # solve for the particle count per particle size
        
        # ----------------------------------------------------------------

        # ----------------------------------------------------------------

        #print(sol_ej.ejecta_ring_bounds[int(ejecta_indices[0]):int(ejecta_indices[np.size(ejecta_indices)-1]), 1])

        # ----------------------------------------------------------------
        # determine the lower and upper bounds of each segment of each of the 20 segments from the data file
        # ----------------------------------------------------------------
        #r_ej_sim_particle_bounds_lower = np.array(sol_ej.ejecta_ring_bounds[int(ejecta_indices[0]):int(ejecta_indices[np.size(ejecta_indices)-1]), 1])
        #r_ej_sim_particle_bounds_upper = np.array(sol_ej.ejecta_ring_bounds[int(ejecta_indices[0]):int(ejecta_indices[np.size(ejecta_indices)-1]), 3])
        lower_bound = int(ejecta_indices[0])
        upper_bound = int(ejecta_indices[np.size(ejecta_indices)-1]+1)
        r_ej_sim_particle_bounds_lower = np.array(sol_ej.ejecta_ring_bounds[lower_bound:upper_bound, 1])
        r_ej_sim_particle_bounds_upper = np.array(sol_ej.ejecta_ring_bounds[lower_bound:upper_bound, 3])
        # ----------------------------------------------------------------
        
        u_ej_sim_particle_mid = np.zeros(np.size(r_ej_sim_particle_bounds_upper)+2)
        r_ej_sim_particle_mid = np.zeros(np.size(r_ej_sim_particle_bounds_upper)+2)
        
        # ----------------------------------------------------------------
        # ----------------------------------------------------------------

        #for i in range(0, np.size(ejecta_velocities[:,1])):
        #    print(i, r_ej_sim_particle_bounds_lower[i], r_ej_sim_particle_bounds_upper[i])



        # ---------------------------------------------------------- 
        # loop through all particle sizes
        # ---------------------------------------------------------- 
        for i in range(0,35):   
        #for i in range(sol_ej.index_grain_sizes[0],sol_ej.index_grain_sizes[1]+1):   

            # solve for the proportion of mass of dust
            proportion = sol_ej.dust_wt_pct[i]

            #print("Particle Size", soil.d_particle[i])
            #print("Proportion", proportion)

            # solve for the volume of the particle
            vol_particle = (4/3) * np.pi * (soil.d_particle[i]/2)**3


            # --------------------------------------------------------------- 
            # extract the ring bounds and the ejecta velocities from datasets
            # --------------------------------------------------------------- 

            # add a beginning and end bounds for the particles ejected for both range and velocity
            # Range: add position 0 and then add the end bound, one more away
            # Velocity: add velocity 0 at both position 0 and then the end bound

            
            for j in range(1, np.size(r_ej_sim_particle_bounds_lower)+1):
                r_ej_sim_particle_mid[j] = (r_ej_sim_particle_bounds_upper[j-1] + r_ej_sim_particle_bounds_lower[j-1])/2
            


            if(int(ejecta_indices[0]) == 0):
                r_ej_sim_particle_mid[0] = 0
            else:
                r_ej_sim_particle_mid[0] =  (sol_ej.ejecta_ring_bounds[int(ejecta_indices[0])-1, 1] + sol_ej.ejecta_ring_bounds[int(ejecta_indices[0]), 1])/2

            r_ej_sim_particle_mid[np.size(r_ej_sim_particle_mid)-1] =  (sol_ej.ejecta_ring_bounds[(int(ejecta_indices[np.size(ejecta_indices)-1])), 3] +  sol_ej.ejecta_ring_bounds[(int(ejecta_indices[np.size(ejecta_indices)-1])+1), 3])/2

            for j in range(1,np.size(u_ej_sim_particle_mid)-1):
                u_ej_sim_particle_mid[j] = ejecta_velocities[j-1,i+1]
            
            #print("mid", r_ej_sim_particle_mid, u_ej_sim_particle_mid)
#

            # interpolate the velocity for the lower and upper bounds
            u_ej_sim_particle_bounds_lower = np.interp(r_ej_sim_particle_bounds_lower, r_ej_sim_particle_mid, u_ej_sim_particle_mid)
            u_ej_sim_particle_bounds_upper = np.interp(r_ej_sim_particle_bounds_upper, r_ej_sim_particle_mid, u_ej_sim_particle_mid)

            #print("Test")
            #print("lower_bound", r_ej_sim_particle_bounds_lower, u_ej_sim_particle_bounds_lower)
            #print("upper_bound", r_ej_sim_particle_bounds_upper, u_ej_sim_particle_bounds_upper)
            

            # ---------------------------------------------------------- 
            # loop through each launch bin range
            # ---------------------------------------------------------- 
            for j in range(0, np.size(u_ej_sim_particle_bounds_lower)):
            
                #print(j, r_ej_sim_particle_bounds_lower[j], r_ej_sim_particle_bounds_upper[j], u_ej_sim_particle_bounds_lower[j], u_ej_sim_particle_bounds_upper[j], ejecta_mass_excavated[j])

                # For each bin, create a large number of simulation particles with a distance from radius and ejecta speed within a bin. 
                r_ej_sim_particle = np.linspace(r_ej_sim_particle_bounds_lower[j], r_ej_sim_particle_bounds_upper[j], sol_ej.n_sim_particles_per_bin)
                u_ej_sim_particle = np.linspace(u_ej_sim_particle_bounds_lower[j], u_ej_sim_particle_bounds_upper[j], sol_ej.n_sim_particles_per_bin)

                #print(r_ej_sim_particle)
                #print(u_ej_sim_particle)

                # compute the total mass and the total number of particles for a certain size excavated from this bin
                total_mass_excavated_particle_size = ejecta_mass_excavated[j] * proportion
                total_count_excavated_particle_size = (total_mass_excavated_particle_size / ejecta_rho_excavated[j]) / vol_particle

                #input("Pause")
                #print(j, ejecta_mass_excavated[j], proportion, total_mass_excavated_particle_size)


                if(r_ej_sim_particle_bounds_lower[j] == 0):
                    lower_bound = 1
                    upper_bound = np.size(r_ej_sim_particle)
                    sim_particles = sol_ej.n_sim_particles_per_bin - 1  
                else:
                    lower_bound = 0
                    upper_bound = np.size(r_ej_sim_particle)
                    sim_particles = sol_ej.n_sim_particles_per_bin
                
                #if(j == 0):
                #    lower_bound = 1
                #    upper_bound = np.size(r_ej_sim_particle)
                #    sim_particles = sol_ej.n_sim_particles_per_bin - 1  

                #elif (j == np.size(u_ej_sim_particle_bounds_lower)-1):
                #    lower_bound = 0
                #    upper_bound = np.size(r_ej_sim_particle) - 1
                #    sim_particles = sol_ej.n_sim_particles_per_bin - 1  

                #else:
                #    lower_bound = 0
                #    upper_bound = np.size(r_ej_sim_particle)
                #    sim_particles = sol_ej.n_sim_particles_per_bin

                # compute the mass of each simulation particle alnong this bin
                mass_sim_particle = total_mass_excavated_particle_size/sim_particles
                count_sim_particle = total_count_excavated_particle_size/sim_particles

               


                #print(" ") 
                #print("Launch Bin #", j)
                #print("Total Mass Excavated Within this Range (kg)", '{:.2e}'.format(ejecta_mass_excavated[j]))
                #print("Mass Excavated for Particle Size ", soil.d_particle[i], "m (kg): ", '{:.2e}'.format(total_mass_excavated_particle_size))
                #print("Count Excavated for Particle Size ", soil.d_particle[i], "m (kg): ", '{:.2e}'.format(total_count_excavated_particle_size))
                #print('Number of Simulation Particles:', sim_particles, ', Mass of a Single Simulation Particle (kg)','{:.2e}'.format(mass_sim_particle))
                #print('Number of Simulation Particles:', sim_particles, ', Count of a Single Simulation Particle ','{:.2e}'.format(count_sim_particle))
                
                #print(" ")

                #time.sleep(1000) 
            
            
                # ---------------------------------------------------------- 
                # loop through each sim particle in each bin
                # ---------------------------------------------------------- 
                #print(u_ej_sim_particle[lower_bound:upper_bound])
                for k in range(lower_bound, upper_bound):

                    escaped_particles = False

                    if(u_ej_sim_particle[k] > 2370):
                        escaped_particles = True
                    else:
                        escaped_particles = False

                        v_ej_bar_2 = u_ej_sim_particle[k]**2 /(1.738E6 * (1.62))
                        r_max = np.abs((2 * 1.738E6 * np.arctan(v_ej_bar_2 * np.sin(3*np.pi/180) * np.cos(3*np.pi/180)/(1 - v_ej_bar_2 * (np.cos(3*np.pi/180))**2))))



                    


                        # determine the maximum arc length range by interpolating the velocities
                        #r_max = np.interp(u_ej_sim_particle[k], range_data[:,0], range_data[:,1])
                        #print(k, u_ej_sim_particle[k], r_max)

                    if(escaped_particles is False):
                        
                        #print(i, j, k, vol_particle, ejecta_rho_excavated[j], u_ej_sim_particle[k], vol_particle * ejecta_rho_excavated[j], 0.5 * vol_particle * ejecta_rho_excavated[j] * (u_ej_sim_particle[k])**2)

                        # sort the sim particle masses into bins based on landing location
                        for l in range(1, sol_ej.n_sorted_bins+1):
                            

                            if(sol_ej.range_bounds[l] > r_max + r_ej_sim_particle[k]):
                                #print(l, sol_ej.range_bounds[l], r_max + r_ej_sim_particle[k])
                                #print(sol_ej.range_bounds[l], r_max + r_ej_sim_particle[k])
                                #print(k, "Belongs in Bin", l - 1, "Particle ", i, r_max, sol_ej.range_bounds[l])

                                if(print_dataset_1):
                                    sol_ej.mass_flux_upon_impact[i,l-1] = sol_ej.mass_flux_upon_impact[i,l-1] + (mass_sim_particle/sol_ej.area_sorted_bins[l-1])
                                    sol_ej.count_flux_upon_impact[i,l-1] = sol_ej.count_flux_upon_impact[i,l-1] + (count_sim_particle/sol_ej.area_sorted_bins[l-1])
                                    sol_ej.mass_upon_impact[i,l-1] = sol_ej.mass_upon_impact[i,l-1] + mass_sim_particle
                                    sol_ej.momentum_flux_upon_impact[i,l-1] = sol_ej.momentum_flux_upon_impact[i,l-1] + (mass_sim_particle * u_ej_sim_particle[k])/sol_ej.area_sorted_bins[l-1]
                                    sol_ej.energy_flux_upon_impact[i,l-1] = sol_ej.energy_flux_upon_impact[i,l-1] + (0.5 * mass_sim_particle * u_ej_sim_particle[k]**2)/sol_ej.area_sorted_bins[l-1]


                                if(print_dataset_2):                                
                                    
                                    energy_particle = 0.5 * vol_particle * ejecta_rho_excavated[j] * (u_ej_sim_particle[k])**2

                                    if(energy_particle >= 1E-7):            
                                    #if(energy_particle >= 1E-9 and energy_particle < 1E-8):            
                                        sol_ej.count_energy_threshold_impact_1[i, l-1] = sol_ej.count_energy_threshold_impact_1[i, l-1] + (count_sim_particle/sol_ej.area_sorted_bins[l-1])
                                    
                                    #if(energy_particle >= 1E-18 and energy_particle < 1E-15):            
                                    if(energy_particle >= 1E-15):            
                                        sol_ej.count_energy_threshold_impact_2[i, l-1] = sol_ej.count_energy_threshold_impact_2[i, l-1] + (count_sim_particle/sol_ej.area_sorted_bins[l-1])
                                    
                                    if(energy_particle >= 1E-12):            
                                        sol_ej.count_energy_threshold_impact_3[i, l-1] = sol_ej.count_energy_threshold_impact_3[i, l-1] + (count_sim_particle/sol_ej.area_sorted_bins[l-1])
                                    
                                    if(energy_particle >= 1E-9):            
                                        sol_ej.count_energy_threshold_impact_4[i, l-1] = sol_ej.count_energy_threshold_impact_4[i, l-1] + (count_sim_particle/sol_ej.area_sorted_bins[l-1])

                                    if(energy_particle >= 1E-6):
                                        sol_ej.count_energy_threshold_impact_5[i, l-1] = sol_ej.count_energy_threshold_impact_5[i, l-1] + (count_sim_particle/sol_ej.area_sorted_bins[l-1])
                                    
                                    if(energy_particle >= 1E-3):
                                        sol_ej.count_energy_threshold_impact_6[i, l-1] = sol_ej.count_energy_threshold_impact_6[i, l-1] + (count_sim_particle/sol_ej.area_sorted_bins[l-1])
                                    
                                    if(energy_particle >= 1E0):
                                        sol_ej.count_energy_threshold_impact_7[i, l-1] = sol_ej.count_energy_threshold_impact_7[i, l-1] + (count_sim_particle/sol_ej.area_sorted_bins[l-1])


                                break
                #input("Pause")
            #print("mass", sol_ej.mass_flux_upon_impact[i,:])
            #print("energy", sol_ej.energy_flux_upon_impact[i,:])
            #print("count", sol_ej.count_flux_upon_impact[i,:])
            #input('press enter to continue')
    if(print_dataset_1):
        np.savetxt(file_folder+identifier+"/"+file_folder+"_post_processing/"+file_folder+"_mass_flux_upon_impact.csv", sol_ej.mass_flux_upon_impact, delimiter=",")
        np.savetxt(file_folder+identifier+"/"+file_folder+"_post_processing/"+file_folder+"_momentum_flux_upon_impact.csv", sol_ej.momentum_flux_upon_impact, delimiter=",")
        np.savetxt(file_folder+identifier+"/"+file_folder+"_post_processing/"+file_folder+"_energy_flux_upon_impact.csv", sol_ej.energy_flux_upon_impact, delimiter=",")
        np.savetxt(file_folder+identifier+"/"+file_folder+"_post_processing/"+file_folder+"_count_flux_upon_impact.csv", sol_ej.count_flux_upon_impact, delimiter=",")
    
    if(print_dataset_2):
        np.savetxt(file_folder+identifier+"/"+file_folder+"_post_processing/"+file_folder+"_count_energy_threshold_impact_1.csv", sol_ej.count_energy_threshold_impact_1, delimiter=",")
        np.savetxt(file_folder+identifier+"/"+file_folder+"_post_processing/"+file_folder+"_count_energy_threshold_impact_2.csv", sol_ej.count_energy_threshold_impact_2, delimiter=",")
        np.savetxt(file_folder+identifier+"/"+file_folder+"_post_processing/"+file_folder+"_count_energy_threshold_impact_3.csv", sol_ej.count_energy_threshold_impact_3, delimiter=",")
        np.savetxt(file_folder+identifier+"/"+file_folder+"_post_processing/"+file_folder+"_count_energy_threshold_impact_4.csv", sol_ej.count_energy_threshold_impact_4, delimiter=",")
        np.savetxt(file_folder+identifier+"/"+file_folder+"_post_processing/"+file_folder+"_count_energy_threshold_impact_5.csv", sol_ej.count_energy_threshold_impact_5, delimiter=",")
        np.savetxt(file_folder+identifier+"/"+file_folder+"_post_processing/"+file_folder+"_count_energy_threshold_impact_6.csv", sol_ej.count_energy_threshold_impact_6, delimiter=",")
        np.savetxt(file_folder+identifier+"/"+file_folder+"_post_processing/"+file_folder+"_count_energy_threshold_impact_7.csv", sol_ej.count_energy_threshold_impact_7, delimiter=",")
#
#plot_deposited_props()

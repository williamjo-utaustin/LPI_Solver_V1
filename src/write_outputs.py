import numpy as np
import matplotlib.pyplot as plt

import var_timestepping as timestep
import var_soil as soil
import sys_output as out
import var_range_of_interest as bounds
import var_impinged_gas as imp
import var_solve_ejection as sol_ej

def print_timestep(ts, h_nozzle):

    print("Writing at time t =", ts * timestep.delta_t, "Nozzle Height: ", h_nozzle)


def write_bounds(bounds_ej_ring):

    if(out.write_bounds):
        np.savetxt("output/ejecta_ring_bounds.csv", bounds_ej_ring, delimiter=',', fmt=' '.join(['%i'] + ['%.4e']*3), header = '# Ejecta Ring Radial Distances Away from Plume Centerline (Index, Min (m), Midpoint (m), Max (m))')
        print("Timestep,", "Depth Excavated,", "Threshold Energy,", "E_down,", "Alpha,","Mdot_flux", "M_area_eroded_inst", "Mdot_cumulative" )

    return None

def plot_erosion_profile(ts, r_from_centerline, d_excavated):
    # compute the total mass eroded in each ring for one timestep
    
    if(out.plot_profile):
    
        if(np.mod(ts,1)==0):
            print("Plotting Figure at Timestep: ", str(ts))
            fig = plt.figure(figsize=(10, 6))
            plt.title("Time = "+str(ts * timestep.delta_t) + " s", fontsize = 18)
            plt.scatter(r_from_centerline, d_excavated)
            #plt.scatter(soil.r_midpoint, soil.h_excavated_mid, linewidth = 2)
            #plt.semilogy(soil.r_midpoint, soil.h_excavated_mid, linewidth = 2)
            plt.xlabel("Distance from Plume Centerline (m)", fontsize = 18)
            plt.ylabel("Excavation Depth (m)", fontsize = 18)
            plt.xticks(fontsize = 18)
            plt.yticks(fontsize = 18)
            #plt.xlim(0,10)
            plt.xlim(0, 100)
            plt.ylim(0.25, 0)
            #plt.gca().set_aspect('equal')
            plt.savefig("output/growth_"+str(ts)+".png", dpi = 300)
            plt.close()

    return None


def write_ejecta_props(ts, h_nozzle, ej_timestep_props, u_ej_timestep, offset_ej_dist_timestep, offset_ej_time_timestep):
    
    # saving values
    if(out.write_props):
        np.savetxt("output/ejecta_properties_"+str(ts)+".csv", ej_timestep_props, delimiter=',', fmt=' '.join(['%i'] + ['%.8e']*4), header = '# header = Nozzle Height (m) ' + str(h_nozzle)+', Offset angle (3 deg), ' + 'Bin Index, Bin Midpoint Location (m), Excavated Density (kg/m^3), Instantaneous Mass Excavated (kg), Total Height Excavated (m)')
        np.savetxt("output/ejecta_velocities_"+str(ts)+".csv", u_ej_timestep, delimiter=',', fmt=' '.join(['%i'] + ['%.4e']*np.size(soil.d_particle)), header = '# Nozzle Height (m) ' + str(h_nozzle)+', Offset angle (3 deg), '+' Ejecta Speeds by Bin Index and Particle Size (O(1), O(10), O(100), O(1000) microns) in sets of 1, 2, 3, 4, 5, 6, 7, 8, 9')
        np.savetxt("output/ejecta_offset_distances_"+str(ts)+".csv", offset_ej_dist_timestep, delimiter=',', fmt=' '.join(['%i'] + ['%.4e']*np.size(soil.d_particle)), header = '# header = Nozzle Height (m) ' + str(h_nozzle)+', Offset angle (3 deg), ' + 'Ejecta Distance Offset by Bin Index and Particle Size (O(1), O(10), O(100), O(1000) microns) in sets of 1, 2, 3, 4, 5, 6, 7, 8, 9')
        np.savetxt("output/ejecta_offset_time_"+str(ts)+".csv", offset_ej_time_timestep, delimiter=',', fmt=' '.join(['%i'] + ['%.4e']*np.size(soil.d_particle)), header = '# header = Nozzle Height (m) ' + str(h_nozzle)+', Offset angle (3 deg), ' + 'Ejecta Time Offset by Bin Index and Particle Size (O(1), O(10), O(100), O(1000) microns) in sets of 1, 2, 3, 4, 5, 6, 7, 8, 9')

    return None

def output_gas_props():
    
    if(out.output_gas):
        for i in range(0,bounds.n_points_centerline-1):
            print(i, soil.r_midpoint[i], imp.v_gas_arr[i], imp.p_gas_arr[i], imp.rho_gas_arr[i], imp.T_gas_arr[i])

    return None

def plot_deposited_props():


    if(sol_ej.plot_only):
        sol_ej.range_bounds_mid = np.genfromtxt("output/range_bounds_"+str(sol_ej.index_grain_sizes[0])+"_"+str(sol_ej.index_grain_sizes[1])+".csv", delimiter=",")
        sol_ej.mass_flux_upon_impact = np.genfromtxt("output/mass_flux_upon_impact_"+str(sol_ej.index_grain_sizes[0])+"_"+str(sol_ej.index_grain_sizes[1])+".csv", delimiter=',')
        sol_ej.momentum_flux_upon_impact = np.genfromtxt("output/momentum_flux_upon_impact_"+str(sol_ej.index_grain_sizes[0])+"_"+str(sol_ej.index_grain_sizes[1])+".csv", delimiter=',')
        sol_ej.energy_flux_upon_impact = np.genfromtxt("output/energy_flux_upon_impact_"+str(sol_ej.index_grain_sizes[0])+"_"+str(sol_ej.index_grain_sizes[1])+".csv", delimiter=',')
        sol_ej.count_flux_upon_impact = np.genfromtxt("output/count_flux_upon_impact_"+str(sol_ej.index_grain_sizes[0])+"_"+str(sol_ej.index_grain_sizes[1])+".csv", delimiter=',')

    total_mass_flux_upon_impact = np.zeros(sol_ej.n_sorted_bins)
    total_momentum_flux_upon_impact = np.zeros(sol_ej.n_sorted_bins)
    total_energy_flux_upon_impact = np.zeros(sol_ej.n_sorted_bins)
    total_count_flux_upon_impact = np.zeros(sol_ej.n_sorted_bins)

    for i in range(0,sol_ej.n_sorted_bins):
        total_mass_flux_upon_impact[i] = np.sum(sol_ej.mass_flux_upon_impact[:,i])
        total_momentum_flux_upon_impact[i] = np.sum(sol_ej.momentum_flux_upon_impact[:,i])
        total_energy_flux_upon_impact[i] = np.sum(sol_ej.energy_flux_upon_impact[:,i])
        total_count_flux_upon_impact[i] = np.sum(sol_ej.count_flux_upon_impact[:,i])

    plt.figure(figsize=(5,5), layout = 'constrained')
    plt.plot(sol_ej.range_bounds_mid, total_energy_flux_upon_impact, color = 'black', linewidth = 2, label = 'Total')

    for i in range(sol_ej.index_grain_sizes[0],sol_ej.index_grain_sizes[1]):
        #plt.ylim(0,10)
        plt.plot(sol_ej.range_bounds_mid, sol_ej.energy_flux_upon_impact[i,:], linewidth = 2, linestyle = 'dotted', label = "Particle Size (m): " + str('{:.1E}'.format(soil.d_particle[i])))
        #plt.xscale('log')
        plt.yscale('log')

    plt.xlabel("Distance from Centerline (m)", fontsize = 16)
    plt.ylabel("Impact Energy $(J/m^2)$", fontsize = 16)
    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)
    plt.locator_params(axis='x', nbins=5)
    plt.xlim(0,sol_ej.range_of_analysis)
    #plt.ylim(1E-1,1E3)
    plt.grid()
    plt.legend(fontsize = 8)
    plt.savefig("output/impact_energy_profile_"+str(sol_ej.index_grain_sizes[0])+"_"+str(sol_ej.index_grain_sizes[1])+".png", bbox_inches='tight', dpi=100)
    plt.close()


    plt.figure(figsize=(5,5), layout = 'constrained')
    plt.plot(sol_ej.range_bounds_mid, total_momentum_flux_upon_impact, color = 'black', linewidth = 2, label = 'Total')

    for i in range(sol_ej.index_grain_sizes[0], sol_ej.index_grain_sizes[1]):
        #plt.ylim(0,10)
        plt.plot(sol_ej.range_bounds_mid, sol_ej.momentum_flux_upon_impact[i,:], linewidth = 2, linestyle = 'dotted', label = "Particle Size (m): " + str('{:.1E}'.format(soil.d_particle[i])))
        #plt.xscale('log')
        plt.yscale('log')

    plt.xlabel("Distance from Centerline (m)", fontsize = 16)
    plt.ylabel("Impact Momentum $(N \cdot s/m^2)$", fontsize = 16)
    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)
    plt.locator_params(axis='x', nbins=5)
    plt.xlim(0,sol_ej.range_of_analysis)
    plt.grid()
    plt.legend(fontsize = 8)
    plt.savefig("output/impact_momentum_profile_"+str(sol_ej.index_grain_sizes[0])+"_"+str(sol_ej.index_grain_sizes[1])+".png", bbox_inches='tight', dpi=100)
    plt.close()

    plt.figure(figsize=(5,5), layout = 'constrained')
    plt.plot(sol_ej.range_bounds_mid, total_mass_flux_upon_impact, color = 'black', linewidth = 2, label = 'Total')

    for i in range(sol_ej.index_grain_sizes[0],sol_ej.index_grain_sizes[1]):
        #plt.ylim(0,10)
        plt.plot(sol_ej.range_bounds_mid, sol_ej.mass_flux_upon_impact[i,:], linewidth = 2, linestyle = 'dotted', label = "Particle Size (m): " + str('{:.1E}'.format(soil.d_particle[i])))
        #plt.xscale('log')
        plt.yscale('log')

    plt.xlabel("Distance from Centerline (m)", fontsize = 16)
    plt.ylabel("Impact Mass $(kg/m^2)$", fontsize = 16)
    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)
    plt.locator_params(axis='x', nbins=5)
    plt.xlim(0,sol_ej.range_of_analysis)
    plt.grid()
    plt.legend(fontsize = 8)
    plt.savefig("output/impact_mass_profile_"+str(sol_ej.index_grain_sizes[0])+"_"+str(sol_ej.index_grain_sizes[1])+".png", bbox_inches='tight',dpi=100)
    plt.close()


    plt.figure(figsize=(5,5), layout = 'constrained')
    plt.plot(sol_ej.range_bounds_mid, total_count_flux_upon_impact, color = 'black', linewidth = 2, label = 'Total')

    for i in range(sol_ej.index_grain_sizes[0],sol_ej.index_grain_sizes[1]):
        #plt.ylim(0,10)
        plt.plot(sol_ej.range_bounds_mid, sol_ej.count_flux_upon_impact[i,:], linewidth = 2, linestyle = 'dotted', label = "Particle Size (m): " + str('{:.1E}'.format(soil.d_particle[i])))
        #plt.xscale('log')
        plt.yscale('log')

    plt.xlabel("Distance from Centerline (m)", fontsize = 16)
    plt.ylabel("Impact Count $(\#/m^2)$", fontsize = 16)
    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)
    plt.locator_params(axis='x', nbins=5)
    plt.xlim(0,sol_ej.range_of_analysis)
    plt.grid()
    plt.legend(fontsize = 8)
    plt.savefig("output/impact_count_profile_"+str(sol_ej.index_grain_sizes[0])+"_"+str(sol_ej.index_grain_sizes[1])+".png", bbox_inches='tight',dpi=100)
    plt.close()


    plt.figure(figsize=(5,5), layout = 'constrained')
    plt.plot(sol_ej.range_bounds_mid, total_count_flux_upon_impact * 0.0025, color = 'black', linewidth = 2, label = 'Total')

    for i in range(sol_ej.index_grain_sizes[0],sol_ej.index_grain_sizes[1]):
        #plt.ylim(0,10)
        plt.plot(sol_ej.range_bounds_mid, sol_ej.count_flux_upon_impact[i,:] * 0.0025, linewidth = 2, linestyle = 'dotted', label = "Particle Size (m): " + str('{:.1E}'.format(soil.d_particle[i])))
        #plt.xscale('log')
        plt.yscale('log')

    plt.xlabel("Distance from Centerline (m)", fontsize = 16)
    plt.ylabel("Impact Count $(\#/25cm^2)$", fontsize = 16)
    plt.xticks(fontsize = 16)
    plt.yticks(fontsize = 16)
    plt.locator_params(axis='x', nbins=5)
    plt.xlim(0,sol_ej.range_of_analysis)
    plt.grid()
    plt.legend(fontsize = 8)
    plt.savefig("output/impact_count_sample_profile_"+str(sol_ej.index_grain_sizes[0])+"_"+str(sol_ej.index_grain_sizes[1])+".png", bbox_inches='tight',dpi=100)
    plt.close()



    return None
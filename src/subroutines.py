import numpy as np
from scipy.integrate import quad
from scipy.optimize import fsolve

import matplotlib.pyplot as plt

import var_range_of_interest as bounds
from fun_impinged_gas import *
from fun_soil import *

def plot_initial_disturbed_altitude():
    
    for j in range(0,bounds.n_points_altitude):

        h_nozzle = bounds.h_nozzle_array[j]
        m_dot = np.zeros(int(bounds.n_points_centerline))

        for i in range(0, bounds.n_points_centerline):

            r_centerline = bounds.r_centerlines_array[i]
            compute_impinged_gas(h_nozzle, r_centerline)
            compute_E_downward() 
            m_dot[i] = compute_mass_erosion_rate()

        plt.plot(bounds.r_centerlines_array, m_dot, linewidth = 3, label = str(h_nozzle) + 'm')

    plt.xlim(0,10000)
    plt.ylim(0,7)
    plt.xticks(fontsize = 18)
    plt.yticks(fontsize = 18)
    plt.xlabel("Distance from Centerline (m)", fontsize = 18)
    plt.ylabel("Unit Area nozzle.Mass Flow Rate ($kg/{m^2s}$)", fontsize = 18)
    plt.legend(title = "Start Altitude (m)", title_fontsize= 16, fontsize = 16 , loc = 'upper right')
    plt.show()


def regolith_removal():

    # initialize_variables()

    # set the timestep (seconds)
    delta_t = 1
    h_nozzle = 4000 #km

    # compute the area of a differential ring
    ring_area = np.zeros(bounds.n_points_centerline-1)
    r_midpoint = np.zeros(bounds.n_points_centerline-1)

    # compute the midpoint distances of each ring
    m_dot_area_eroded = np.zeros(bounds.n_points_centerline)
    m_dot_eroded = np.zeros(bounds.n_points_centerline-1)
    m_excavated_inst = np.zeros(bounds.n_points_centerline-1)
    m_excavated_cumulative = np.zeros(bounds.n_points_centerline-1)
    
    # initialize excavation heights 
    h_excavated = np.zeros(bounds.n_points_centerline-1)
    h_excavated_old = np.zeros(bounds.n_points_centerline-1)
    surf_den_excavated = np.zeros(bounds.n_points_centerline-1)
    
    for i in range(1,bounds.n_points_centerline):
        r_midpoint[i-1] = (bounds.r_centerlines_array[i] + bounds.r_centerlines_array[i-1])/2
    
    # compute the ring areas
    for i in range(1,bounds.n_points_centerline):
        ring_area[i-1] = np.pi * ((bounds.r_centerlines_array[i])**2 - (bounds.r_centerlines_array[i-1])**2) 
        #print(i-1, "to", i, bounds.r_centerlines_array[i-1], "to", bounds.r_centerlines_array[i], ring_area[i-1])


    # loop()

    # --------------------------------------
    # begin looping over here (ie timesteps)
    # --------------------------------------

    # change nozzle height (v_descent & delta_t)

    # solve for the mass erosion rate by area (kg/m^2 s) 
    for i in range(0, bounds.n_points_centerline):
        r_centerline = bounds.r_centerlines_array[i]
        compute_impinged_gas(h_nozzle, r_centerline)
        compute_E_downward() 
        m_dot_area_eroded[i] = compute_mass_erosion_rate()
    
    # compute the ring mass erosion rate (kg/s) at each timestep
    for i in range(1,bounds.n_points_centerline):
        m_dot_eroded[i-1] = ((m_dot_area_eroded[i] + m_dot_area_eroded[i-1])/2) * ring_area[i-1]
    
    # compute the total mass eroded in each ring for one timestep
    m_excavated_inst = m_dot_eroded * delta_t

    # compute the cumulative total mass eroded in each ring
    m_excavated_cumulative = m_excavated_cumulative + m_excavated_inst

    # compute for the depth of the new regolith
    for i in range(0,bounds.n_points_centerline-1):
        
        # save excavation depth for current timestep
        h_excavated_old[i] = h_excavated[i]

        # solve for the surficial density (kg/m^2)
        surf_den_excavated[i] = (m_excavated_inst[i])/ring_area[i]
        
        # solve for the depth z of excavation by integration and solving for the upper bound
        # the lower bound is defined as the starting depth at each timestep 
        def func(z):
            y, err = quad(soil_density_depth, h_excavated[i], z)
            return surf_den_excavated[i] - y
        sol = fsolve(func, 1E-8)

        # update the new excavated depth
        h_excavated[i] = sol[0]

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
    #plt.semilogy(r_midpoint, h_excavated)
    #plt.xlim(bounds.min_centerline, bounds.max_centerline)
    #plt.ylim(1E-8, 1E1)
    #plt.show()
        
    return None
import numpy as np
from scipy.integrate import quad
from scipy.optimize import fsolve

import matplotlib.pyplot as plt

import var_timestepping as timestep
import var_range_of_interest as bounds
import var_soil as soil
from fun_impinged_gas import *
from fun_soil import *
from fun_nozzle import *

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


def regolith_removal():

    # initialize_variables()
    # start from t = 0 or from a restart file

    # unload the restart file


    # set the timestep (seconds)
    h_nozzle = 4000 #km

    # --------------------------------------
    # begin looping over here (ie timesteps)
    # --------------------------------------

    print("Timestep,", "Depth Excavated,", "Threshold Energy,", "E_down,", "Alpha,","Mdot")
    for t in range(0,5000):
        
        compute_mass_eroded_quantities(h_nozzle)

        # need to turn this into a function but it doesn't seem to be working
        for i in range(0,bounds.n_points_centerline-1):

            # save excavation depth for current timestep
            soil.h_excavated_mid_old[i] = soil.h_excavated_mid[i]

            # solve for the surficial density (kg/m^2)
            soil.surf_den_excavated[i] = (soil.m_excavated_inst[i])/soil.ring_area[i]

            # solve for the depth z of excavation by integration and solving for the upper bound
            def func(z):
                y, err = quad(compute_soil_density_depth, soil.h_excavated_mid[i], z)
                return soil.surf_den_excavated[i] - y
            
            sol = fsolve(func, 1E-8)

            soil.h_excavated_mid[i] = sol[0]

        set_new_excavation_depths()

        # update the new nozzle velocity
        h_nozzle = update_nozzle_height(h_nozzle)
        
        # print statisics
        if(np.mod(t,10)==0):
            print(t, soil.h_excavated_bounds[860], soil.E_th[860], imp.E_downward, soil.alpha[860], soil.m_dot_area_eroded[860])

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
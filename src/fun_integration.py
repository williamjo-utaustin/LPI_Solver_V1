import numpy as np
import var_constants as cs
import var_range_of_interest as bounds
import var_soil as soil

import var_impinged_gas as imp

from fun_soil import *
from fun_integration import *



def cd_sphere(Re, Kn):
    
    if Re <= 0.0:
        CD = 0.0
    elif Re > 8.0e6:
        CD = 0.2
    elif Re > 0.0 and Re <= 2:
        CD = (24.0/Re)
    elif Re > 2 and Re <= 100.0:
        p = np.array([4.22, -14.05, 34.87, 0.658])
        CD = np.polyval(p, 1.0/Re) 
    elif Re > 100.0 and Re <= 1.0e4:
        p = np.array([-30.41, 43.72, -17.08, 2.41])
        CD = np.polyval(p, 1.0/np.log10(Re))
    elif Re > 1.0e4 and Re <= 3.35e5:
        p = np.array([-0.1584, 2.031, -8.472, 11.932])
        CD = np.polyval(p, np.log10(Re))
    elif Re > 3.35e5 and Re <= 5.0e5:
        x1 = np.log10(Re/4.5e5)
        CD = 91.08*x1**4 + 0.0764
    else:
        p = np.array([-0.06338, 1.1905, -7.332, 14.93])
        CD = np.polyval(p, np.log10(Re))

    # add in correction factor to account for rarefied effects
    S_correction = 1 + Kn * (2.514 + 0.8 * np.exp(-0.55/Kn))
    CD = CD/S_correction
    return CD

def dudt(u_p, t, d_p, rho_p, u_g, p_g, rho_g, T_g):
    
    A_const= 1.71575E-7
    beta = 0.78
    mfp = np.sqrt(np.pi/(2 * p_g * rho_g)) * A_const * T_g**beta
    
    Kn = mfp/d_p
    Re = (rho_g * np.abs(u_g - u_p) * d_p)/(A_const * T_g ** beta)
    
    if Re <= 0.0:
        CD = 0.0
    elif Re > 0.0 and Re <= 2:
        CD = (24.0/Re)
    elif Re >= 2 and Re < 500:
        CD = 18.5 * Re **-0.6
    else:
        CD = 0.44

    # add in correction factor to account for rarefied effects
    S_correction = 1 + Kn * (2.514 + 0.8 * np.exp(-0.55/Kn))
    CD = CD/S_correction

    C_constant = 3 * rho_g * CD / (4 * rho_p * d_p)
    return C_constant * (u_g - u_p)**2


def determine_ejecta_props(t):
    
    # copy the contents of density array (will be the previous timestep)
    avg_rho_p_excavated_old = soil.avg_rho_p_excavated.copy()

    # create an index count, and counter for total affected indices
    index_count = np.zeros(bounds.n_points_centerline-1)
    total_affected_indices = 0

    # compute the average density of the material that has already been excavated
    # loop through all bins in the centerline array
    for i in range(0,bounds.n_points_centerline-1):
        
        # if the excavated height at this current timestep is deeper than the previous timestep (has been dug)
        # compute the average density of the soil from the depth dug
        if((soil.h_excavated_mid[i] - soil.h_excavated_mid_old[i]) > 0):
            soil.avg_rho_p_excavated[i],err = quad(compute_soil_density_depth, soil.h_excavated_mid_old[i], soil.h_excavated_mid[i]) * 1/(soil.h_excavated_mid[i] - soil.h_excavated_mid_old[i])
        
        # compute the density at the excavated height
        else:
            soil.avg_rho_p_excavated[i] = compute_soil_density_depth(soil.h_excavated_mid[i])  

        # if the bin has been excavated at this timestep, there should be a change in height
        if(np.abs(soil.avg_rho_p_excavated[i] - avg_rho_p_excavated_old[i] > 0)):
            
            # write down the bin number in a list called 'index count'
            # move to the next index in index count
            index_count[total_affected_indices] = i
            total_affected_indices = total_affected_indices + 1
   
    rho_gas_arr_mid = np.zeros(bounds.n_points_centerline-1)
    u_gas_arr_mid = np.zeros(bounds.n_points_centerline-1)
    T_gas_arr_mid = np.zeros(bounds.n_points_centerline-1)
    P_gas_arr_mid = np.zeros(bounds.n_points_centerline-1)


    u_ej_timestep = np.zeros([total_affected_indices-1, np.size(soil.d_particle) + 1])
    offset_ej_dist_timestep = np.zeros([total_affected_indices-1, np.size(soil.d_particle) + 1])
    offset_ej_time_timestep = np.zeros([total_affected_indices-1, np.size(soil. d_particle) + 1])
    ej_timestep_props = np.zeros([total_affected_indices-1, 5])

    for i in range(0,total_affected_indices-1):

        i_n = int(index_count[i])
        ip1 = int(index_count[i+1])

        ej_timestep_props[i, 0] = i_n
        ej_timestep_props[i, 1] = soil.r_midpoint[i_n]
        ej_timestep_props[i, 2] = soil.avg_rho_p_excavated[i_n]
        ej_timestep_props[i, 3] = soil.m_excavated_inst[i_n]
        ej_timestep_props[i ,4] = soil.h_excavated_mid[i_n]


        print("Integrating from index ", i_n, "to ", ip1, "Iteration #", i+1,"out of ",total_affected_indices-1)

        rho_gas_arr_mid[i_n] = (imp.rho_gas_arr[i_n] + imp.rho_gas_arr[ip1])/2
        u_gas_arr_mid[i_n] = (imp.v_gas_arr[i_n] + imp.v_gas_arr[ip1])/2
        T_gas_arr_mid[i_n] = (imp.T_gas_arr[i_n] + imp.T_gas_arr[ip1])/2
        P_gas_arr_mid[i_n] = (imp.p_gas_arr[i_n] + imp.p_gas_arr[ip1])/2


        #print(u_gas_arr_mid[i_n])

        u_ej_timestep[i,0] = i_n
        offset_ej_dist_timestep[i,0] = i_n
        offset_ej_time_timestep[i,0] = i_n

        for j in range(0,np.size(soil.d_particle)):
            
            d_p = soil.d_particle[j]
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
                    offset_ej_time_timestep[i,j+1] = t_sub
                    break

                delta_t = du/(dudt(u_p, ts, d_p, soil.avg_rho_p_excavated[i_n], u_gas_arr_mid[i_n], P_gas_arr_mid[i_n], rho_gas_arr_mid[i_n], T_gas_arr_mid[i_n]))

                u_p = u_p + du
                t_sub = t_sub + delta_t

                x_p = x_p + (u_p * delta_t) * np.cos(3 * np.pi/180)
                y_p = y_p + (u_p * delta_t) * np.sin(3 * np.pi/180)


    return u_ej_timestep, offset_ej_dist_timestep, offset_ej_time_timestep, ej_timestep_props

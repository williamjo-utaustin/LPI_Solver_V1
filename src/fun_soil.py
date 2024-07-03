import numpy as np

from scipy.integrate import quad
from scipy.optimize import fsolve

import var_range_of_interest as bounds

import var_soil as soil
import var_impinged_gas as imp
import var_constants as cs
import var_timestepping as timestep

from fun_impinged_gas import *

def compute_soil_density_depth(z):
    #return soil.rho_inf * (z + 0.122) / (z + 0.18)
    soil.rho_bulk = soil.rho_0 - (soil.rho_0 - soil.rho_inf)*(1-np.exp(-z/soil.h)) 
    return soil.rho_bulk

def compute_cohesive_energy_density(z):
    return soil.alpha_0 * np.exp((compute_soil_density_depth(z) - soil.rho_0)/soil.k_cohesion)

def compute_threshold_energy(z):
    E_th = soil.E_th_0 * compute_ratio(z)
    return E_th

def compute_E_downward():
    imp.E_downward = imp.rho_gas * imp.v_gas**2 * imp.v_T /12
    return None

def compute_ratio(z):    
    ratio = (soil.rho_0 * cs.g * soil.D_bracket + soil.alpha_0)/(compute_soil_density_depth(z) * cs.g * soil.D_bracket + compute_cohesive_energy_density(z))
    return ratio


def compute_mass_erosion_rate(z):
    
    E_th = compute_threshold_energy(z)
    alpha = compute_cohesive_energy_density(z)

    compute_E_downward()
    
    if(imp.E_downward < E_th):
        m_dot = 0
    else:
        m_dot = soil.rho_0 * soil.epsilon * (imp.E_downward - soil.E_th_0)/(soil.rho_0 * cs.g * 1.5 * soil.D_bracket + soil.alpha_0)

    m_dot = m_dot * compute_ratio(z)

    return m_dot, E_th, alpha


def compute_mass_eroded_quantities(h_nozzle):
    
    n_bins_eroded_dt = 0
    bin_eroded_id = np.zeros(bounds.n_points_centerline)

    for i in range(0, bounds.n_points_centerline):
        
        initial_m_dot_area_eroded = soil.m_dot_area_eroded[i]

        r_centerline = bounds.r_centerlines_array[i]
        compute_impinged_gas(h_nozzle, r_centerline)

        imp.rho_gas_arr[i] = imp.rho_gas        
        imp.p_gas_arr[i] = imp.p_gas
        imp.v_gas_arr[i] = imp.v_gas      
        imp.T_gas_arr[i] = imp.T_gas        
        imp.v_T_arr[i] = imp.v_T        

        soil.m_dot_area_eroded[i], soil.E_th[i], soil.alpha[i] = compute_mass_erosion_rate(soil.h_excavated_bounds[i])


        print(i, r_centerline, imp.v_gas_arr[i], soil.m_dot_area_eroded[i], imp.E_downward, soil.E_th[i], soil.alpha[i])

        if (soil.m_dot_area_eroded[i] - initial_m_dot_area_eroded > 0):
            
            bin_eroded_id[n_bins_eroded_dt] = i
            n_bins_eroded_dt = n_bins_eroded_dt + 1

    # we only need to loop these variables instead of the whole loop 
    #print(n_bins_eroded_dt, bin_eroded_id)
    
    # compute the ring mass erosion rate (kg/s) at each timestep at the midpoints
    #for i in range(1,n_bins_eroded_dt):
    #    
    #    i_m1 = int(bin_eroded_id[i-1])
    #    i_p0 = int(bin_eroded_id[i])

    #    soil.m_dot_eroded_mid[i_m1] = ((soil.m_dot_area_eroded[i_p0] + soil.m_dot_area_eroded[i_m1])/2) * soil.ring_area[i_m1]

    for i in range(0, bounds.n_points_centerline-1):
        soil.m_dot_eroded_mid[i] = ((soil.m_dot_area_eroded[i] + soil.m_dot_area_eroded[i+1])/2) * soil.ring_area[i]
        #print("test", i, soil.m_dot_eroded_mid[i], soil.m_dot_eroded_mid[i]/soil.ring_area[i])


    # compute the total mass eroded in each ring for one timestep
    soil.m_excavated_inst = soil.m_dot_eroded_mid * timestep.delta_t

    # compute the cumulative total mass eroded in each ring
    soil.m_excavated_cumulative = soil.m_excavated_cumulative + soil.m_excavated_inst
    
    #for i in range(0, bounds.n_points_centerline-1):
    #    print("test", i, soil.m_dot_eroded_mid[i], soil.m_dot_eroded_mid[i]/soil.ring_area[i], soil.m_excavated_inst[i], soil.m_excavated_cumulative[i])

    return None

def set_new_excavation_depths():
    soil.h_excavated_bounds[0] = soil.h_excavated_mid[0]
    for i in range(1, bounds.n_points_centerline-1):
        soil.h_excavated_bounds[i] = (soil.h_excavated_mid[i-1] + soil.h_excavated_mid[i])/2
    soil.h_excavated_bounds[bounds.n_points_centerline-1] = soil.h_excavated_mid[bounds.n_points_centerline-2]
    
    return None


def integral_equation(h_new_exc, h_old_exc, mass_excavated):
    integral_equation = mass_excavated - quad(compute_soil_density_depth, h_old_exc, h_new_exc)[0]
    return integral_equation


def solve_excavation_depths():
        
    for i in range(0,bounds.n_points_centerline-1):

        # save excavation depth for current timestep
        soil.h_excavated_mid_old[i] = soil.h_excavated_mid[i]

        # solve for the surficial density (kg/m^2)
        soil.surf_den_excavated[i] = (soil.m_excavated_inst[i])/soil.ring_area[i]
        soil.h_excavated_mid[i] = fsolve(integral_equation, soil.h_excavated_mid_old[i], args = (soil.h_excavated_mid_old[i], soil.surf_den_excavated[i]))

    return None
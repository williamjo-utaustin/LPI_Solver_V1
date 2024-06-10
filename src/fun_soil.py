import numpy as np

from scipy.integrate import quad
from scipy.optimize import fsolve

import var_range_of_interest as bounds
import var_excavation as excavate

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
    
    for i in range(0, bounds.n_points_centerline):
        r_centerline = bounds.r_centerlines_array[i]
        compute_impinged_gas(h_nozzle, r_centerline)
        soil.m_dot_area_eroded[i], soil.E_th[i], soil.alpha[i] = compute_mass_erosion_rate(soil.h_excavated_bounds[i])
    
    # compute the ring mass erosion rate (kg/s) at each timestep at the midpoints
    for i in range(1,bounds.n_points_centerline):
        soil.m_dot_eroded_mid[i-1] = ((soil.m_dot_area_eroded[i] + soil.m_dot_area_eroded[i-1])/2) * soil.ring_area[i-1]

    # compute the total mass eroded in each ring for one timestep
    soil.m_excavated_inst = soil.m_dot_eroded_mid * timestep.delta_t

    # compute the cumulative total mass eroded in each ring
    soil.m_excavated_cumulative = soil.m_excavated_cumulative + soil.m_excavated_inst

    return None

def set_new_excavation_depths():
    soil.h_excavated_bounds[0] = soil.h_excavated_mid[0]
    for i in range(1, bounds.n_points_centerline-1):
        soil.h_excavated_bounds[i] = (soil.h_excavated_mid[i-1] + soil.h_excavated_mid[i])/2
    soil.h_excavated_bounds[bounds.n_points_centerline-1] = soil.h_excavated_mid[bounds.n_points_centerline-2]
    return None


import numpy as np
import var_soil as soil
import var_impinged_gas as imp
import var_constants as cs

def compute_soil_density_depth(z):
    #return soil.rho_inf * (z + 0.122) / (z + 0.18)
    soil.rho_bulk = soil.rho_0 - (soil.rho_0 - soil.rho_inf)*(1-np.exp(-z/soil.h)) 
    return soil.rho_bulk

def compute_cohesive_energy_density(z):
    return soil.alpha_0 * np.exp((compute_soil_density_depth(z) - soil.rho_inf)/soil.k_cohesion)

def compute_threshold_energy(z):
    E_th = soil.E_th_0 * (soil.rho_0 * cs.g * soil.D_84 + soil.alpha_0)/(compute_soil_density_depth(z)*cs.g*soil.D_84 + compute_cohesive_energy_density(z))
    return E_th

def compute_E_downward():
    imp.E_downward = imp.rho_gas * imp.v_gas**2 * imp.v_T /12
    return None

def compute_mass_erosion_rate(E_th, alpha):
    if(imp.E_downward < E_th):
        m_dot = 0
    else:
        m_dot = soil.rho_bulk * soil.epsilon * (imp.E_downward - E_th)/(soil.rho_bulk * cs.g * 1.5 * soil.D_84 + alpha)
    return m_dot
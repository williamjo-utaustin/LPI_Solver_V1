import numpy as np
import var_soil as soil
import var_impinged_gas as imp
import var_constants as cs

def soil_density_depth(z):
    return soil.rho_inf * (z + 0.122) / (z + 0.18)

def cohesive_energy_density(z):
    alpha_0 = 0.289
    rho_0 = soil.rho_inf * (0.122/0.18)
    return alpha_0 * np.exp((soil_density_depth(z) - soil.rho_inf)/soil.k_cohesion)

def compute_E_downward():
    imp.E_downward = imp.rho_gas * imp.v_gas**2 * imp.v_T /12
    return None

def compute_mass_erosion_rate():
    if(imp.E_downward < soil.E_th):
        m_dot = 0
    else:
        m_dot = soil.rho_bulk_soil * soil.epsilon * (imp.E_downward - soil.E_th)/(soil.rho_bulk_soil * cs.g * 1.5 * soil.D_84 + soil.alpha)
    return m_dot
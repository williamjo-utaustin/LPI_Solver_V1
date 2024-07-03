import numpy as np
import var_nozzle as nozzle
import var_impinged_gas as imp
import var_constants as cs
from fun_compressible import pNS_over_p0

def compute_impinged_gas(h_nozzle, r_centerline):

    # solve for theta
    theta = np.arctan(r_centerline/h_nozzle)
    imp.rho_gas = nozzle.rho_e * (nozzle.k_bar/2) * (h_nozzle/(nozzle.r_nozzle))**(-2) * (np.cos(theta))**nozzle.k_bar
    p_normal_shock = nozzle.P_0 * pNS_over_p0(nozzle.Ma,nozzle.gamma)
    
    # compute the gas stagnation pressure at the surface if the nozzle is far enough away
    surface_criteria = (h_nozzle/(nozzle.r_nozzle)) * np.sqrt(2/(nozzle.k_bar+2))

    if(surface_criteria > 1):
        p_stag = p_normal_shock * ((nozzle.k_bar+2)/2) * (h_nozzle/(nozzle.r_nozzle))**-2
    else:
        p_stag = p_normal_shock

    xi = np.sqrt((nozzle.k_bar + 4)/2) * np.tan(theta)
    # compute the gas pressure 
    #imp.p_gas = p_stag * np.exp(-(r_centerline/h_nozzle)**2)
    imp.p_gas = p_stag * np.exp(-(xi)**2)
    #imp.p_gas = nozzle.throttle /((np.pi * h_nozzle**2) * (2/(nozzle.k_bar + 2)))

    # compute the gas velocity (From Morris Dissertation)
    imp.v_gas = np.sqrt((2/imp.rho_gas) * p_stag * (xi)**2 * np.exp(-(xi)**2))
    #imp.v_gas = np.sqrt(((2 * nozzle.gamma * cs.g * nozzle.R_gas * nozzle.T_0)/(nozzle.gamma - 1)) * (1 - np.cos(theta) ** ((nozzle.gamma-1)*(nozzle.k_bar+4)/nozzle.gamma)))
    #imp.v_gas = np.sqrt(((2 * nozzle.gamma * cs.g * nozzle.R_gas * nozzle.T_0)/(nozzle.gamma - 1)))
    #imp.v_gas = np.sqrt(2 * (xi)**2 * nozzle.R_gas * nozzle.T_0)
 
    #imp.v_gas = np.sqrt((2/imp.rho_gas) * p_stag * (r_centerline/h_nozzle)**2 * ((nozzle.k_bar+4)/2) * np.exp(-((nozzle.k_bar+4)/2) * (r_centerline/h_nozzle)**2))
    #imp.v_gas = np.sqrt(2 * (r_centerline/h_nozzle)**2 * nozzle.R_gas * nozzle.T_0)

    ## this value in Arnowitz is supposed to be a power (See Eq 57 and Eq 60)
    imp.T_gas = nozzle.T_0 * np.cos(theta) ** ((nozzle.gamma-1)*(nozzle.k_bar+4)/nozzle.gamma)
    #imp.T_gas = nozzle.T_0 * np.cos(theta) 
    #imp.T_gas = nozzle.T_0
    #imp.T_gas = nozzle.T_0 * np.cos(theta) * ((nozzle.gamma-1)/nozzle.gamma)
    #imp.T_gas = imp.p_gas/(imp.rho_gas * nozzle.R_gas) 
    #print(T_gas) 

    #imp.v_T = np.sqrt((2/imp.rho_gas) * p_stag * (xi)**2 * np.exp(-(xi)**2))
    #imp.v_T = np.sqrt(2 * (xi)**2 * nozzle.R_gas * nozzle.T_0)
    #imp.v_T = np.sqrt(2 * (r_centerline/h_nozzle)**2 * nozzle.R_gas * nozzle.T_0)
    imp.v_T = np.sqrt(nozzle.R_gas * imp.T_gas)

    print("--------------------------------------------------------------------------------------------------------------")
    print("For a height", str(h_nozzle), "meters and distance away from nozzle centerline", str(r_centerline), "meters")
    print("Xi: ", xi)
    print("Gas Pressure (Pa): ", imp.p_gas)
    print("Gas Temperature (K): ", imp.T_gas)
    print("Gas Density (kg/m^3): ", imp.rho_gas)
    print("Gas Velocity (m/s): ", imp.v_gas)
    print("Gas Mean Thermal Velocity (m/s): ", imp.v_T)
    print("--------------------------------------------------------------------------------------------------------------")

    return None
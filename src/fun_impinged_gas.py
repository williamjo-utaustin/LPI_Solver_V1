import numpy as np
import var_nozzle as nozzle
import var_impinged_gas as imp

from fun_compressible import pNS_over_p0

def compute_impinged_gas(h_nozzle, r_centerline):

    theta = np.arctan(r_centerline/h_nozzle)
    imp.rho_gas = nozzle.rho_e * (nozzle.k_bar/2) * (h_nozzle/(nozzle.D_nozzle/2))**(-2) * (np.cos(theta))**nozzle.k_bar
    p_normal_shock = nozzle.P_0 * pNS_over_p0(nozzle.Ma,nozzle.gamma)
        # compute the gas stagnation pressure at the surface if the nozzle is far enough away
    surface_criteria = (h_nozzle/(nozzle.D_nozzle/2)) * np.sqrt(2/(nozzle.k_bar+2))

    if(surface_criteria > 1):
        p_stag = p_normal_shock * ((nozzle.k_bar+2)/2) * (h_nozzle/(nozzle.D_nozzle/2))**-2
    else:
        p_stag = p_normal_shock

    # compute the gas pressure 
    imp.p_gas = p_stag * np.exp(-(r_centerline/h_nozzle)**2)
    #print(nozzle.P_0, p_nornozzle.Mal_shock, p_stag, p_gas)

    # compute the gas velocity
    imp.v_gas = np.sqrt((2/imp.rho_gas) * p_stag * (r_centerline/h_nozzle)**2 * ((nozzle.k_bar+4)/2) * np.exp(-((nozzle.k_bar+4)/2) * (r_centerline/h_nozzle)**2))
    #print(v_gas)

    ## this value is good
    imp.T_gas = nozzle.T_0 * np.cos(theta) * (nozzle.gamma-1)*(nozzle.k_bar+4)/nozzle.gamma
    #print(T_gas) 

    imp.v_T = np.sqrt(2 * (r_centerline/h_nozzle)**2 * nozzle.R_gas * nozzle.T_0)

    #print("--------------------------------------------------------------------------------------------------------------")
    #print("For a height", str(h_nozzle), "meters and distance away from nozzle centerline", str(r_centerline), "meters")
    #print("Gas Pressure (Pa): ", p_gas)
    #print("Gas Temperature (K): ", T_gas)
    #print("Gas Density (kg/m^3): ", rho_gas)
    #print("Gas Velocity (m/s): ", v_gas)
    #print("Gas Mean Thernozzle.Mal Velocity (m/s): ", v_T)
    #print("--------------------------------------------------------------------------------------------------------------")

    return None
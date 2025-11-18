import numpy as np

import sys
#sys.path.append('/Users/williamjo/Documents/LPI/Codes/plume_regolith_solver_v1/src')

import var_constants as cs
import var_nozzle as nozzle
import var_timestepping as timestep
import spacecraft_options as spacecraft
from fun_compressible import *

def compute_gas_Mw():
    # insert dataset for gas molar mass
    gas_data = np.genfromtxt('data/FontesEtAl_GasData.csv',delimiter = ',')

    # compute molar mass from gas data
    M_w = 0
    for i in range(0,np.size(gas_data[:,0])):
        M_w = M_w + gas_data[i,0] * gas_data[i,1]

    return M_w/1000 # convert to kg/mol

def compute_nozzle_exhaust():
    #nozzle.m_dot_e = mDot_over_A(nozzle.Ma, nozzle.gamma, nozzle.P_0, nozzle.R_gas, nozzle.T_0) * nozzle.A_nozzle
    #nozzle.P_e = nozzle.P_0/(p0_over_p(nozzle.Ma, nozzle.gamma))
    #nozzle.T_e = nozzle.T_0/(T0_over_T(nozzle.Ma, nozzle.gamma))
    #nozzle.rho_e = (nozzle.P_0/(nozzle.R_gas*nozzle.T_0)) / (rho0_over_rho(nozzle.Ma, nozzle.gamma))
    #nozzle.v_e = nozzle.Ma * np.sqrt(nozzle.gamma * nozzle.R_gas * nozzle.T_e)

    # Use the CEA Calculator from NASA GRC to Find Exhaust Properties
    if(spacecraft.lander_type == 'starship_raptor_nominal_50' or spacecraft.lander_type == 'starship_raptor_vel_high_50'):
        
        # new values (December 2024)
        nozzle.m_dot_e = 14.519 * 3 #6.9556 * 3
        nozzle.P_e = 273.05 #1.9201E-3 * 100000 
        nozzle.T_e = 895.39 #1934.6
        nozzle.rho_e = 6.9627E-4 #2.2507E-4
        nozzle.v_e = nozzle.Ma * np.sqrt(nozzle.gamma * nozzle.R_gas * nozzle.T_e)
       
        # old values (July 2024)
        #nozzle.m_dot_e = 6.9556 * 3 #6.9556 * 3
        #nozzle.P_e = 150.81 #1.9201E-3 * 100000 
        #nozzle.T_e = 1129.3 #1934.6
        #nozzle.rho_e = 3.049E-4 #2.2507E-4
        #nozzle.v_e = nozzle.Ma * np.sqrt(nozzle.gamma * nozzle.R_gas * nozzle.T_e)

    elif(spacecraft.lander_type == 'starship_raptor_nominal_100' or spacecraft.lander_type == 'starship_raptor_vel_high_100' or spacecraft.lander_type == 'starship_raptor_vel_mid_100'  or spacecraft.lander_type == 'starship_raptor_vel_low_100'):
        
        # new values (December 2024)
        nozzle.m_dot_e = 28.977 * 3 #6.9556 * 3
        nozzle.P_e = 542.8 #1.9201E-3 * 100000 
        nozzle.T_e = 890.24 #1934.6
        nozzle.rho_e = 1.3921E-2 #2.2507E-4
        nozzle.v_e = nozzle.Ma * np.sqrt(nozzle.gamma * nozzle.R_gas * nozzle.T_e)
       
        # old values (July 2024)
        #nozzle.m_dot_e = 13.879 * 3 #12.121958 * 3
        #nozzle.P_e = 297.43 #3.7318E-3 * 100000 
        #nozzle.T_e = 1113.2 #1904.5
        #nozzle.rho_e = 6.100E-4 #4.4565E-4
        #nozzle.v_e = nozzle.Ma * np.sqrt(nozzle.gamma * nozzle.R_gas * nozzle.T_e)
    
    elif(spacecraft.lander_type == 'starship_thruster_nominal_50' or spacecraft.lander_type == 'starship_thruster_vel_high_50'):
        nozzle.m_dot_e = 2.9510 * 9 #2.69027987885064 * 9
        nozzle.P_e = 1400 #1.4168E-02 * 100000 
        nozzle.T_e = 2291 #2560.7
        nozzle.rho_e = 1.387E-3 #1.1203E-3
        nozzle.v_e = nozzle.Ma * np.sqrt(nozzle.gamma * nozzle.R_gas * nozzle.T_e)

    elif(spacecraft.lander_type == 'starship_thruster_nominal_100' or spacecraft.lander_type == 'starship_thruster_vel_high_100' or spacecraft.lander_type == 'starship_thruster_vel_mid_100'  or spacecraft.lander_type == 'starship_thruster_vel_low_100'):
        nozzle.m_dot_e = 5.855 * 9 #5.325825686994666 * 9
        nozzle.P_e = 2777.9 #2.8198E-2 * 100000 
        nozzle.T_e = 2317.6 #2622.8
        nozzle.rho_e = 2.667E-3 #2.1961E-3
        nozzle.v_e = nozzle.Ma * np.sqrt(nozzle.gamma * nozzle.R_gas * nozzle.T_e)

    elif('apollo_' in spacecraft.lander_type):
        nozzle.m_dot_e = 3.3809
        nozzle.P_e = 209.23
        nozzle.T_e = 1042.9
        nozzle.rho_e = 5.1237E-4
        nozzle.v_e = nozzle.Ma * np.sqrt(nozzle.gamma * nozzle.R_gas * nozzle.T_e)

    else:
        exit("No Specified Lander Type")

    #if(spacecraft.lander_type == 'starship_low_thrusters_vel_nom'):
    #    nozzle.m_dot_e = 7.5722372611511615
    #    nozzle.P_e = 2.3672E-3 * 100000 
    #    nozzle.T_e = 1926.3
    #    nozzle.rho_e = 2.7895E-4
    #    nozzle.v_e = 4.2663 * np.sqrt(1.1696 * 440.5393 * nozzle.T_e)
    #if(spacecraft.lander_type == 'starship_low_thrusters_vel_high'):
    #    nozzle.m_dot_e = 7.5722372611511615
    #    nozzle.P_e = 2.3672E-3 * 100000 
    #    nozzle.T_e = 1926.3
    #    nozzle.rho_e = 2.7895E-4
    #    nozzle.v_e = 4.2663 * np.sqrt(1.1696 * 440.5393 * nozzle.T_e)
    #if(spacecraft.lander_type == 'starship_low_thrusters_vel_med'):
    #    nozzle.m_dot_e = 7.5722372611511615
    #    nozzle.P_e = 2.3672E-3 * 100000 
    #    nozzle.T_e = 1926.3
    #    nozzle.rho_e = 2.7895E-4
    #    nozzle.v_e = 4.2663 * np.sqrt(1.1696 * 440.5393 * nozzle.T_e)
    #if(spacecraft.lander_type == 'starship_low_thrusters_vel_low'):
    #    nozzle.m_dot_e = 7.5722372611511615
    #    nozzle.P_e = 2.3672E-3 * 100000 
    #    nozzle.T_e = 1926.3
    #    nozzle.rho_e = 2.7895E-4
    #    nozzle.v_e = 4.2663 * np.sqrt(1.1696 * 440.5393 * nozzle.T_e)

    return None

# WJ new notes
# 10/11/25
# set an upward accelleration to counteract the negative velocity
# need to make sure that the height will be 0 and the will be 0, even with displacement
# IE maybe we need to make a variable that just tracks from the base and add displacement for the nozzle exhaust height
# basically we need to replace 1.5 in nozzle.v_descent and in h_nozzle with a user inputted accelleration (a_thrust)
# Need to also figure out how to limit the timesteps in the code (ie not all things go past 100 s in descent time)

def update_nozzle_height(h_nozzle, timestep_count):
   
    print(spacecraft.lander_type)

    if(spacecraft.lander_type == 'starship_thruster_nominal_100' or spacecraft.lander_type == 'starship_thruster_nominal_50'):
        if(h_nozzle > 233):
            nozzle.v_descent = 11
        elif(h_nozzle <= 233 and h_nozzle > 183):
            nozzle.v_descent = 8.5
        elif(h_nozzle <= 183 and h_nozzle > 133):
            nozzle.v_descent = 6
        elif(h_nozzle <= 133 and h_nozzle > 83):
            nozzle.v_descent = 2.5
        elif(h_nozzle <= 83 and h_nozzle > 33):
            nozzle.v_descent = 1
        elif(h_nozzle <= 33 and h_nozzle > 0):
            nozzle.v_descent = 0.7
        else:
            nozzle.v_descent = 0.7
    
    elif(spacecraft.lander_type == 'starship_raptor_nominal_100' or spacecraft.lander_type == 'starship_raptor_nominal_50'):
        
        
        # new model (December 2024) 
        
        nozzle.v_descent = nozzle.v_init - 1.5 * (timestep.delta_t * timestep_count)
        h_nozzle = nozzle.h_nozzle_init - nozzle.v_init * (timestep.delta_t * timestep_count) + 0.5 * 1.5 * (timestep.delta_t * timestep_count)**2

        # old model (July 2024)
        #if(h_nozzle > 251):
        #    nozzle.v_descent = 13
        #elif(h_nozzle <= 251 and h_nozzle > 201):
        #    nozzle.v_descent = 11
        #elif(h_nozzle <= 201 and h_nozzle > 151):
        #    nozzle.v_descent = 8.5
        #elif(h_nozzle <= 151 and h_nozzle > 101):
        #    nozzle.v_descent = 6
        #elif(h_nozzle <= 101 and h_nozzle > 51):
        #    nozzle.v_descent = 2.5
        #elif(h_nozzle <= 51 and h_nozzle > 1):
        #    nozzle.v_descent = 1
        #else:
        #    nozzle.v_descent = 0.7
    
    elif('apollo_' in spacecraft.lander_type):
        pass    
    
    elif('starship_thruster_vel_' in spacecraft.lander_type):
        pass
    elif('starship_raptor_vel_' in spacecraft.lander_type):
        pass

    else:
        exit("Not a Specified Lander Type")
    
    
    # update the new nozzle velocity
    # new model (December 2024)
    h_nozzle = h_nozzle + 1

    # old model (July 2024)
    #h_nozzle = h_nozzle - (nozzle.v_descent * timestep.delta_t)
    return h_nozzle

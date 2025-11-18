import numpy as np
from fun_nozzle import compute_gas_Mw
import var_constants as cs
from fun_compressible import *
import spacecraft_options as spacecraft

# ---------------------------------------------------------
# Change only these variables (Input Conditions)
# ---------------------------------------------------------
#gamma = 1.2
#P_0 = 19046 # Pascals (11 Ton)
##P_0 = 69259.33 # Pascals (40 Ton)
#T_0 = 4875 # Kelvin
#Ma = 2.5
#D_nozzle = 1.6 # meters
# ------------------------------------------------------------------------

a_thrust = 1.5




# ---------------------------------------------------------
# new parameters (12/9/2024)
# ---------------------------------------------------------
if(spacecraft.lander_type == 'starship_raptor_nominal_50'):
    P_0 =  475000 
    T_0 = 3012
    rho_0 = 0.34811
    gamma = 1.2086
    Ma = 5.055
    D_nozzle = 4.745819212  # 2.74 m diameter for 1 engine, we have 3
    v_descent = 14
    v_init = v_descent
    h_nozzle_init = 65
    min_altitude_lander = 1
    R_gas = 437.97
    

if(spacecraft.lander_type == 'starship_raptor_nominal_100'):
    P_0 = 955000 
    T_0 = 3069 
    rho_0 = 0.6905 
    gamma = 1.208 
    Ma = 5.0759 
    D_nozzle = 4.745819212  # 2.74 m diameter for 1 engine, we have 3
    v_descent = 14
    v_init = v_descent
    h_nozzle_init = 65
    min_altitude_lander = 1
    R_gas = 437.97 

# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------

#if(spacecraft.lander_type == 'starship_raptor_nominal_50'):
#    P_0 =  235570 #216560 
#    T_0 = 3140 #3423.5
#    rho_0 = 0.15848 #0.11988
#    gamma = 1.2231 #1.1636
#    Ma = 4.8265 #4.255
#    D_nozzle = 4.745819212  # 2.74 m diameter for 1 engine, we have 3
#    v_descent = 11
#    h_nozzle_init = 251
#    min_altitude_lander = 1
#    R_gas = 438.00355 #440.9761010244944

if(spacecraft.lander_type == 'starship_raptor_vel_high_50'):
    P_0 =  235570 #216560 
    T_0 = 3140 #3423.5
    rho_0 = 0.15848 #0.11988
    gamma = 1.2231 #1.1636
    Ma = 4.8265 #4.255
    D_nozzle = 4.745819212  # 2.74 m diameter for 1 engine, we have 3
    v_descent = 11
    h_nozzle_init = 251
    min_altitude_lander = 1
    R_gas = 438.00355 #440.9761010244944

#if(spacecraft.lander_type == 'starship_raptor_nominal_100'):
#    P_0 = 476390 #436920.0 
#    T_0 = 3225.2 #3542.8
#    rho_0 = 0.31447 #0.23609
#    gamma = 1.2224 #1.1824
#    Ma = 4.8762 #4.2948
#    D_nozzle = 4.745819212  # 2.74 m diameter for 1 engine, we have 3
#    v_descent = 11
#    h_nozzle_init = 251
#    min_altitude_lander = 1
#    R_gas = 437.9646 #439.6868453626443

if(spacecraft.lander_type == 'starship_raptor_vel_high_100'):
    P_0 = 476390 #436920.0 
    T_0 = 3225.2 #3542.8
    rho_0 = 0.31447 #0.23609
    gamma = 1.2224 #1.1824
    Ma = 4.8762 #4.2948
    D_nozzle = 4.745819212  # 2.74 m diameter for 1 engine, we have 3
    v_descent = 11
    h_nozzle_init = 251
    min_altitude_lander = 1
    R_gas = 437.9646 #439.6868453626443

if(spacecraft.lander_type == 'starship_raptor_vel_mid_100'):
    P_0 = 476390 #436920.0 
    T_0 = 3225.2 #3542.8
    rho_0 = 0.31447 #0.23609
    gamma = 1.2224 #1.1824
    Ma = 4.8762 #4.2948
    D_nozzle = 4.745819212  # 2.74 m diameter for 1 engine, we have 3
    v_descent = 6
    h_nozzle_init = 251
    min_altitude_lander = 1
    R_gas = 437.9646 #439.6868453626443

if(spacecraft.lander_type == 'starship_raptor_vel_low_100'):
    P_0 = 476390 #436920.0 
    T_0 = 3225.2 #3542.8
    rho_0 = 0.31447 #0.23609
    gamma = 1.2224 #1.1824
    Ma = 4.8762 #4.2948
    D_nozzle = 4.745819212  # 2.74 m diameter for 1 engine, we have 3
    v_descent = 1
    h_nozzle_init = 251
    min_altitude_lander = 1
    R_gas = 437.9646 #439.6868453626443

if(spacecraft.lander_type == 'starship_thruster_nominal_50'):
    P_0 = 27259 #26986.0
    T_0 = 2886 #3101.3 
    rho_0 = 1.94E-2 #1.605E-2
    gamma = 1.1243 #1.09999
    Ma = 2.4868 #2.4858
    D_nozzle = 3  # 1 m diameter for 1 engine, we have 9 engines
    v_descent = 11
    h_nozzle_init = 283
    min_altitude_lander = 33
    R_gas = 449.886 #493.873

if(spacecraft.lander_type == 'starship_thruster_vel_high_50'):
    P_0 = 27259 #26986.0
    T_0 = 2886 #3101.3 
    rho_0 = 1.94E-2 #1.605E-2
    gamma = 1.1243 #1.09999
    Ma = 2.4868 #2.4858
    D_nozzle = 3  # 1 m diameter for 1 engine, we have 9 engines
    v_descent = 11
    h_nozzle_init = 283
    min_altitude_lander = 33
    R_gas = 449.886 #493.873


if(spacecraft.lander_type == 'starship_thruster_nominal_100'):
    P_0 = 54689 #54059.0
    T_0 = 2966 #3203.6 
    rho_0 = 0.038336 #0.031394
    gamma = 1.1327 #1.1031
    Ma = 2.4922 #2.4905
    D_nozzle = 3  # 2.74 m diameter for 1 engine, we have 3
    v_descent = 11 
    h_nozzle_init = 283
    min_altitude_lander = 33
    R_gas = 447.643 #489.55

if(spacecraft.lander_type == 'starship_thruster_vel_high_100'):
    P_0 = 54689 #54059.0
    T_0 = 2966 #3203.6 
    rho_0 = 0.038336 #0.031394
    gamma = 1.1327 #1.1031
    Ma = 2.4922 #2.4905
    D_nozzle = 3  # 2.74 m diameter for 1 engine, we have 3
    v_descent = 11 
    h_nozzle_init = 283
    min_altitude_lander = 33
    R_gas = 447.643 #489.55

if(spacecraft.lander_type == 'starship_thruster_vel_mid_100'):
    P_0 = 54689 #54059.0
    T_0 = 2966 #3203.6 
    rho_0 = 0.038336 #0.031394
    gamma = 1.1327 #1.1031
    Ma = 2.4922 #2.4905
    D_nozzle = 3  # 2.74 m diameter for 1 engine, we have 3
    v_descent = 6 
    h_nozzle_init = 283
    min_altitude_lander = 33
    R_gas = 447.643 #489.55
if(spacecraft.lander_type == 'starship_thruster_vel_low_100'):
    P_0 = 54689 #54059.0
    T_0 = 2966 #3203.6 
    rho_0 = 0.038336 #0.031394
    gamma = 1.1327 #1.1031
    Ma = 2.4922 #2.4905
    D_nozzle = 3  # 2.74 m diameter for 1 engine, we have 3
    v_descent = 1
    h_nozzle_init = 283
    min_altitude_lander = 33
    R_gas = 447.643 #489.55

#if(spacecraft.lander_type == "apollo_vel_slow_alt_high"):
if('apollo_' in spacecraft.lander_type):
    P_0 = 154320 #54059.0
    T_0 = 2956 #3203.6 
    rho_0 = 0.12656 #0.031394
    gamma = 1.2575 #1.1031
    Ma = 4.507 #2.4905
    D_nozzle = 1.6  # 2.74 m diameter for 1 engine, we have 3

    min_altitude_lander = 1
    R_gas = 391.559 #489.55
    
    if(spacecraft.lander_type=="apollo_vel_fast_alt_mid"):
        print("writing conditions for apollo_vel_fast_alt_mid")
        v_descent = 11 
        h_nozzle_init = 32.5
    elif(spacecraft.lander_type=="apollo_vel_slow_alt_mid"):
        print("writing conditions for apollo_vel_slow_alt_mid")
        v_descent = 1 
        h_nozzle_init = 32.5
    elif(spacecraft.lander_type=="apollo_vel_mid_alt_mid"):
        print("writing conditions for apollo_vel_mid_alt_mid")
        v_descent = 6
        h_nozzle_init = 32.5
    elif(spacecraft.lander_type=="apollo_vel_slow_alt_high"):
        print("writing conditions for apollo_vel_slow_alt_high")
        v_descent = 1 
        h_nozzle_init = 45
    elif(spacecraft.lander_type=="apollo_vel_slow_alt_low"):
        print("writing conditions for apollo_vel_slow_alt_low")
        v_descent = 1
        h_nozzle_init = 15
    else:
        exit("Apollo descent speed not found!")














#if(spacecraft.lander_type == 'apollo'):
#    gamma = 1.2
#    P_0 = 179000
#    T_0 = 1612
#    Ma = 4.359
#    D_nozzle = 1.5 
#    v_descent = 0.7
#    h_nozzle_init = 31.5
#    min_altitude_lander = 1
#    M_w = compute_gas_Mw()
#    
#    # determine gas/nozzle variables at the exit
#    R_gas = cs.R_universal/M_w # J/(kgK)
#    
#elif(spacecraft.lander_type == 'starship_nominal'):
#    gamma = 1.2
#    P_0 = 69259.33
#    T_0 = 4875
#    Ma = 2.5
#    D_nozzle = 1.6 
#    r_nozzle = D_nozzle/2
#    v_descent = 11
#    h_nozzle_init = 283 # starship nozzle height plus 5 times lander height
#    min_altitude_lander = 33 # maximum nozzle height
#    M_w = compute_gas_Mw()
#    
#    # determine gas/nozzle variables at the exit
#    R_gas = cs.R_universal/M_w # J/(kgK)
#elif(spacecraft.lander_type == 'starship_vel_high'):
#    gamma = 1.2
#    P_0 = 69259.33
#    T_0 = 4875
#    Ma = 2.5
#    D_nozzle = 1.6 
#    r_nozzle = D_nozzle/2
#    v_descent = 11
#    h_nozzle_init = 283 # starship nozzle height plus 5 times lander height
#    min_altitude_lander = 33 # maximum nozzle height
#    M_w = compute_gas_Mw()
#    
#    # determine gas/nozzle variables at the exit
#    R_gas = cs.R_universal/M_w # J/(kgK)
#elif(spacecraft.lander_type == 'starship_vel_low'):
#    gamma = 1.2
#    P_0 = 69259.33
#    T_0 = 4875
#    Ma = 2.5
#    D_nozzle = 1.6 
#    r_nozzle = D_nozzle/2
#    v_descent = 1
#    h_nozzle_init = 283 # starship nozzle height plus 5 times lander height
#    min_altitude_lander = 33 # maximum nozzle height
#    M_w = compute_gas_Mw()
#    
#    # determine gas/nozzle variables at the exit
#    R_gas = cs.R_universal/M_w # J/(kgK)
#elif(spacecraft.lander_type == 'starship_vel_mid'):
#    gamma = 1.2
#    P_0 = 69259.33
#    T_0 = 4875
#    Ma = 2.5
#    D_nozzle = 1.6 
#    r_nozzle = D_nozzle/2
#    v_descent = 6
#    h_nozzle_init = 283 # starship nozzle height plus 5 times lander height
#    min_altitude_lander = 33 # maximum nozzle height
#    
#    M_w = compute_gas_Mw()
#    
#    # determine gas/nozzle variables at the exit
#    R_gas = cs.R_universal/M_w # J/(kgK)
#
#elif(spacecraft.lander_type == 'starship_mid_thrusters_vel_nom'):
#    gamma = 1.09
#    P_0 = 33500
#    T_0 = 3122.56
#    Ma = 2.481
#    D_nozzle = 3.0 
#    v_descent = 11
#    h_nozzle_init = 283 # starship nozzle height plus 5 times lander height
#    min_altitude_lander = 33 # maximum nozzle height
#    
#    # determine gas/nozzle variables at the exit
#    R_gas = 420.829 # J/(kgK)
#elif(spacecraft.lander_type == 'starship_mid_thrusters_vel_high'):
#    gamma = 1.09
#    P_0 = 33500
#    T_0 = 3122.56
#    Ma = 2.481
#    D_nozzle = 3.0 
#    v_descent = 11
#    h_nozzle_init = 283 # starship nozzle height plus 5 times lander height
#    min_altitude_lander = 33 # maximum nozzle height
#    # determine gas/nozzle variables at the exit
#    R_gas = 420.829 # J/(kgK)
#elif(spacecraft.lander_type == 'starship_mid_thrusters_vel_low'):
#    gamma = 1.09
#    P_0 = 33500
#    T_0 = 3122.56
#    Ma = 2.481
#    D_nozzle = 3.0 
#    v_descent = 1
#    h_nozzle_init = 283 # starship nozzle height plus 5 times lander height
#    min_altitude_lander = 33 # maximum nozzle height
#    # determine gas/nozzle variables at the exit
#    R_gas = 420.829 # J/(kgK)
#elif(spacecraft.lander_type == 'starship_mid_thrusters_vel_mid'):
#    gamma = 1.09
#    P_0 = 33500
#    T_0 = 3122.56
#    Ma = 2.481
#    D_nozzle = 3.0 
#    v_descent = 6
#    h_nozzle_init = 283 # starship nozzle height plus 5 times lander height
#    min_altitude_lander = 33 # maximum nozzle height
#    # determine gas/nozzle variables at the exit
#    R_gas = 420.829 # J/(kgK)
#
#elif(spacecraft.lander_type == 'starship_low_thrusters_vel_nom'):
#    gamma = 1.1696
#    P_0 = 270000
#    T_0 = 3460
#    Ma = 4.2663
#    D_nozzle = 4.7458 
#    v_descent = 11
#    h_nozzle_init = 250 # starship nozzle height plus 5 times lander height
#    min_altitude_lander = 1 # maximum nozzle height
#    
#    # determine gas/nozzle variables at the exit
#    R_gas = 440.5393 # J/(kgK)
#elif(spacecraft.lander_type == 'starship_low_thrusters_vel_high'):
#    gamma = 1.1696
#    P_0 = 270000
#    T_0 = 3460
#    Ma = 4.2663
#    D_nozzle = 3.0 
#    v_descent = 11
#    h_nozzle_init = 250 # starship nozzle height plus 5 times lander height
#    min_altitude_lander = 1 # maximum nozzle height
#    # determine gas/nozzle variables at the exit
#    R_gas = 440.5393 # J/(kgK)
#elif(spacecraft.lander_type == 'starship_low_thrusters_vel_low'):
#    gamma = 1.1696
#    P_0 = 270000
#    T_0 = 3460
#    Ma = 4.2663
#    D_nozzle = 3.0 
#    v_descent = 1
#    h_nozzle_init = 250 # starship nozzle height plus 5 times lander height
#    min_altitude_lander = 1 # maximum nozzle height
#    # determine gas/nozzle variables at the exit
#    R_gas = 440.5393 # J/(kgK)
#elif(spacecraft.lander_type == 'starship_low_thrusters_vel_mid'):
#    gamma = 1.1696
#    P_0 = 270000
#    T_0 = 3460
#    Ma = 4.2663
#    D_nozzle = 3.0 
#    v_descent = 6
#    h_nozzle_init = 250 # starship nozzle height plus 5 times lander height
#    min_altitude_lander = 1 # maximum nozzle height
#    # determine gas/nozzle variables at the exit
#    R_gas = 440.5393 # J/(kgK)
#
#else:
#    print("Error! Need to specify a valid spacecraft type.")
#    exit()

r_nozzle = D_nozzle/2


 # calculate nozzle properties
A_nozzle = (np.pi/4) * D_nozzle**2 #m^2

if("starship_raptor_" in spacecraft.lander_type):
    A_throat = A_nozzle/107
elif("starship_thruster_" in spacecraft.lander_type):
    A_throat = A_nozzle/4
elif("apollo_" in spacecraft.lander_type):
    A_throat = A_nozzle/53.6
else:
    exit("Incorrect Lander Name Leader")
#A_throat = A_nozzle/A_over_Astar(Ma, gamma) # m^2
#
#if(spacecraft.lander_type == 'starship_low_thrusters_vel_nom'):
#    A_throat = A_nozzle / 107
#if(spacecraft.lander_type == 'starship_low_thrusters_vel_high'):
#    A_throat = A_nozzle / 107
#if(spacecraft.lander_type == 'starship_low_thrusters_vel_mid'):
#    A_throat = A_nozzle / 107
#if(spacecraft.lander_type == 'starship_low_thrusters_vel_low'):
#    A_throat = A_nozzle / 107
#if(spacecraft.lander_type == 'starship_metholox_nominal_100' or spacecraft.lander_type == 'starship_metholox_nominal_50'):
#    A_throat = A_nozzle / 4



k_bar = gamma * (gamma - 1) * Ma**2

# initialize exit conditions
m_dot_e = None
P_e = None 
T_e = None 
rho_e = None
v_e = None

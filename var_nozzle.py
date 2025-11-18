import numpy as np
from fun_nozzle import compute_gas_Mw
import var_constants as cs
from fun_compressible import *
import spacecraft_options as spacecraft
import var_timestepping as timesteps

# Declare variables that would be used later
v_descent = None

if spacecraft.lander_type == 'bluemoon_H30_lambda1p01_m1p5e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 15000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 0.694 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 30 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 60.8581 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 0.9859006 # m/s (downward)
    a_thrust = 0.0162 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  458552.92 # Pa
    T_0 = 3204.7876 # K
    rho_0 = 0.24059298 # kg/m3

    # List Exit Conditions
    P_e = 535.19526 # Pa
    T_e = 1543.6204 # K
    rho_e = 0.00064743508 # kg/m^3
    m_dot_e = 1.8307779 * n_engines #kg/s
    gamma = 1.2218855
    Ma = 4.2446737

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'bluemoon_H30_lambda1p5_m1p5e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 15000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 0.694 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 30 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 8.60663 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 6.97137 # m/s (downward)
    a_thrust = 0.81 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  682350.16 # Pa
    T_0 = 3256.3053 # K
    rho_0 = 0.35422271 # kg/m3

    # List Exit Conditions
    P_e = 781.94578 # Pa
    T_e = 1526.5296 # K
    rho_e = 0.00095631794 # kg/m^3
    m_dot_e = 2.7130224 * n_engines #kg/s
    gamma = 1.2235941
    Ma = 4.2774114

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'bluemoon_H30_lambda2_m1p5e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 15000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 0.694 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 30 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 6.08581 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 9.859006 # m/s (downward)
    a_thrust = 1.62 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  911141.54 # Pa
    T_0 = 3293.8264 # K
    rho_0 = 0.46937191 # kg/m3

    # List Exit Conditions
    P_e = 1030.7217 # Pa
    T_e = 1514.7321 # K
    rho_e = 0.0012702528 # kg/m^3
    m_dot_e = 3.612218 * n_engines #kg/s
    gamma = 1.2246446
    Ma = 4.3003136

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'bluemoon_H30_lambda2p5_m1p5e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 15000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 0.694 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 30 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 4.96904 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 12.074767 # m/s (downward)
    a_thrust = 2.43 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  1140249.8 # Pa
    T_0 = 3322.9482 # K
    rho_0 = 0.58393856 # kg/m3

    # List Exit Conditions
    P_e = 1277.2478 # Pa
    T_e = 1505.9718 # K
    rho_e = 0.0015834008 # kg/m^3
    m_dot_e = 4.5109027 * n_engines #kg/s
    gamma = 1.2254207
    Ma = 4.3175967

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'bluemoon_H30_lambda3_m1p5e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 15000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 0.694 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 30 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 4.30331 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 13.94274 # m/s (downward)
    a_thrust = 3.24 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  1369629.5 # Pa
    T_0 = 3346.7844 # K
    rho_0 = 0.69806498 # kg/m3

    # List Exit Conditions
    P_e = 1522.2554 # Pa
    T_e = 1498.9504 # K
    rho_e = 0.0018958658 # kg/m^3
    m_dot_e = 5.4086981 * n_engines #kg/s
    gamma = 1.2259393
    Ma = 4.3314993

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'bluemoon_H30_lambda4_m1p5e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 15000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 0.694 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 30 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 3.51364 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 17.076299 # m/s (downward)
    a_thrust = 4.86 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  1829007.5 # Pa
    T_0 = 3384.5595 # K
    rho_0 = 0.92528781 # kg/m3

    # List Exit Conditions
    P_e = 2008.5449 # Pa
    T_e = 1488.3198 # K
    rho_e = 0.0025193762 # kg/m^3
    m_dot_e = 7.2029886 * n_engines #kg/s
    gamma = 1.226758
    Ma = 4.3527313

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'bluemoon_H30_lambda5_m1p5e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 15000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 0.694 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 30 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 3.0429 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 19.718012 # m/s (downward)
    a_thrust = 6.48 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  2289051.9 # Pa
    T_0 = 3413.865 # K
    rho_0 = 1.1514354 # kg/m3

    # List Exit Conditions
    P_e = 2491.1483 # Pa
    T_e = 1480.3752 # K
    rho_e = 0.0031414201 # kg/m^3
    m_dot_e = 8.9954512 * n_engines #kg/s
    gamma = 1.2273891
    Ma = 4.3687606

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'bluemoon_H60_lambda1p5_m1p5e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 15000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 0.694 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 60 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 12.1716 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 9.859006 # m/s (downward)
    a_thrust = 0.81 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  682350.16 # Pa
    T_0 = 3256.3053 # K
    rho_0 = 0.35422271 # kg/m3

    # List Exit Conditions
    P_e = 781.94578 # Pa
    T_e = 1526.5296 # K
    rho_e = 0.00095631794 # kg/m^3
    m_dot_e = 2.7130224 * n_engines #kg/s
    gamma = 1.2235941
    Ma = 4.2774114

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'bluemoon_H60_lambda2_m1p5e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 15000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 0.694 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 60 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 8.60663 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 13.94274 # m/s (downward)
    a_thrust = 1.62 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  911141.54 # Pa
    T_0 = 3293.8264 # K
    rho_0 = 0.46937191 # kg/m3

    # List Exit Conditions
    P_e = 1030.7217 # Pa
    T_e = 1514.7321 # K
    rho_e = 0.0012702528 # kg/m^3
    m_dot_e = 3.612218 * n_engines #kg/s
    gamma = 1.2246446
    Ma = 4.3003136

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'bluemoon_H60_lambda2p5_m1p5e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 15000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 0.694 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 60 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 7.02728 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 17.076299 # m/s (downward)
    a_thrust = 2.43 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  1140249.8 # Pa
    T_0 = 3322.9482 # K
    rho_0 = 0.58393856 # kg/m3

    # List Exit Conditions
    P_e = 1277.2478 # Pa
    T_e = 1505.9718 # K
    rho_e = 0.0015834008 # kg/m^3
    m_dot_e = 4.5109027 * n_engines #kg/s
    gamma = 1.2254207
    Ma = 4.3175967

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'bluemoon_H60_lambda3_m1p5e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 15000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 0.694 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 60 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 6.08581 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 19.718012 # m/s (downward)
    a_thrust = 3.24 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  1369629.5 # Pa
    T_0 = 3346.7844 # K
    rho_0 = 0.69806498 # kg/m3

    # List Exit Conditions
    P_e = 1522.2554 # Pa
    T_e = 1498.9504 # K
    rho_e = 0.0018958658 # kg/m^3
    m_dot_e = 5.4086981 * n_engines #kg/s
    gamma = 1.2259393
    Ma = 4.3314993

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'bluemoon_H60_lambda4_m1p5e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 15000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 0.694 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 60 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 4.96904 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 24.149534 # m/s (downward)
    a_thrust = 4.86 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  1829007.5 # Pa
    T_0 = 3384.5595 # K
    rho_0 = 0.92528781 # kg/m3

    # List Exit Conditions
    P_e = 2008.5449 # Pa
    T_e = 1488.3198 # K
    rho_e = 0.0025193762 # kg/m^3
    m_dot_e = 7.2029886 * n_engines #kg/s
    gamma = 1.226758
    Ma = 4.3527313

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'bluemoon_H60_lambda5_m1p5e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 15000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 0.694 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 60 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 4.30331 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 27.88548 # m/s (downward)
    a_thrust = 6.48 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  2289051.9 # Pa
    T_0 = 3413.865 # K
    rho_0 = 1.1514354 # kg/m3

    # List Exit Conditions
    P_e = 2491.1483 # Pa
    T_e = 1480.3752 # K
    rho_e = 0.0031414201 # kg/m^3
    m_dot_e = 8.9954512 * n_engines #kg/s
    gamma = 1.2273891
    Ma = 4.3687606

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'bluemoon_H30_lambda1p5_m3e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 30000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 0.694 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 30 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 8.60663 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 6.97137 # m/s (downward)
    a_thrust = 0.81 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  1369629.5 # Pa
    T_0 = 3346.7844 # K
    rho_0 = 0.69806498 # kg/m3

    # List Exit Conditions
    P_e = 1522.2554 # Pa
    T_e = 1498.9504 # K
    rho_e = 0.0018958658 # kg/m^3
    m_dot_e = 5.4086981 * n_engines #kg/s
    gamma = 1.2259393
    Ma = 4.3314993

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'bluemoon_H30_lambda2_m3e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 30000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 0.694 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 30 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 6.08581 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 9.859006 # m/s (downward)
    a_thrust = 1.62 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  1829007.5 # Pa
    T_0 = 3384.5595 # K
    rho_0 = 0.92528781 # kg/m3

    # List Exit Conditions
    P_e = 2008.5449 # Pa
    T_e = 1488.3198 # K
    rho_e = 0.0025193762 # kg/m^3
    m_dot_e = 7.2029886 * n_engines #kg/s
    gamma = 1.226758
    Ma = 4.3527313

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'bluemoon_H30_lambda2p5_m3e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 30000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 0.694 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 30 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 4.96904 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 12.074767 # m/s (downward)
    a_thrust = 2.43 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  2289051.9 # Pa
    T_0 = 3413.865 # K
    rho_0 = 1.1514354 # kg/m3

    # List Exit Conditions
    P_e = 2491.1483 # Pa
    T_e = 1480.3752 # K
    rho_e = 0.0031414201 # kg/m^3
    m_dot_e = 8.9954512 * n_engines #kg/s
    gamma = 1.2273891
    Ma = 4.3687606

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'bluemoon_H60_lambda1p5_m3e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 30000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 0.694 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 60 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 12.1716 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 9.859006 # m/s (downward)
    a_thrust = 0.81 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  1369629.5 # Pa
    T_0 = 3346.7844 # K
    rho_0 = 0.69806498 # kg/m3

    # List Exit Conditions
    P_e = 1522.2554 # Pa
    T_e = 1498.9504 # K
    rho_e = 0.0018958658 # kg/m^3
    m_dot_e = 5.4086981 * n_engines #kg/s
    gamma = 1.2259393
    Ma = 4.3314993

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'bluemoon_H60_lambda2_m3e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 30000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 0.694 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 60 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 8.60663 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 13.94274 # m/s (downward)
    a_thrust = 1.62 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  1829007.5 # Pa
    T_0 = 3384.5595 # K
    rho_0 = 0.92528781 # kg/m3

    # List Exit Conditions
    P_e = 2008.5449 # Pa
    T_e = 1488.3198 # K
    rho_e = 0.0025193762 # kg/m^3
    m_dot_e = 7.2029886 * n_engines #kg/s
    gamma = 1.226758
    Ma = 4.3527313

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'bluemoon_H60_lambda2p5_m3e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 30000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 0.694 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 60 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 7.02728 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 17.076299 # m/s (downward)
    a_thrust = 2.43 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  2289051.9 # Pa
    T_0 = 3413.865 # K
    rho_0 = 1.1514354 # kg/m3

    # List Exit Conditions
    P_e = 2491.1483 # Pa
    T_e = 1480.3752 # K
    rho_e = 0.0031414201 # kg/m^3
    m_dot_e = 8.9954512 * n_engines #kg/s
    gamma = 1.2273891
    Ma = 4.3687606

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'bluemoon_H30_lambda1p5_m4p5e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 45000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 0.694 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 30 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 8.60663 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 6.97137 # m/s (downward)
    a_thrust = 0.81 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  2058928 # Pa
    T_0 = 3400.0125 # K
    rho_0 = 1.0384926 # kg/m3

    # List Exit Conditions
    P_e = 2250.2155 # Pa
    T_e = 1484.0982 # K
    rho_e = 0.0028305296 # kg/m^3
    m_dot_e = 8.0993142 * n_engines #kg/s
    gamma = 1.2271179
    Ma = 4.3612625

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'bluemoon_H60_lambda1p5_m4p5e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 45000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 0.694 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 60 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 12.1716 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 9.859006 # m/s (downward)
    a_thrust = 0.81 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  2058928 # Pa
    T_0 = 3400.0125 # K
    rho_0 = 1.0384926 # kg/m3

    # List Exit Conditions
    P_e = 2250.2155 # Pa
    T_e = 1484.0982 # K
    rho_e = 0.0028305296 # kg/m^3
    m_dot_e = 8.0993142 * n_engines #kg/s
    gamma = 1.2271179
    Ma = 4.3612625

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

# --------------------------------------------------------------------------
# ADDED NOV 1st, 2025
# --------------------------------------------------------------------------
if spacecraft.lander_type == 'bluemoon_H30_lambda1p5_m2p25e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 22500 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 0.694 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 30 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 8.60663 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 6.97137 # m/s (downward)
    a_thrust = 0.81 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  1025641.4 # Pa
    T_0 = 3309.1795 # K
    rho_0 = 0.52670559 # kg/m3

    # List Exit Conditions
    P_e = 1154.1646 # Pa
    T_e = 1510.0513 # K
    rho_e = 0.0014269262 # kg/m^3
    m_dot_e = 4.061754 * n_engines #kg/s
    gamma = 1.2250769
    Ma = 4.3094718

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'bluemoon_H30_lambda1p5_m3p75e04':     
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 37500 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 0.694 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 30 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 8.60663 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 6.97137 # m/s (downward)
    a_thrust = 0.81 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  1714083 # Pa
    T_0 = 3376.0562 # K
    rho_0 = 0.86858998 # kg/m3

    # List Exit Conditions
    P_e = 1887.3576 # Pa
    T_e = 1490.7043 # K
    rho_e = 0.0023636825 # kg/m^3
    m_dot_e = 6.7544524 * n_engines #kg/s
    gamma = 1.2266141
    Ma = 4.3480774

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)



if spacecraft.lander_type == 'bluemoon_H60_lambda1p5_m2p25e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 22500 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 0.694 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 60 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 12.1716 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 9.859006 # m/s (downward)
    a_thrust = 0.81 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  1025641.4 # Pa
    T_0 = 3309.1795 # K
    rho_0 = 0.52670559 # kg/m3

    # List Exit Conditions
    P_e = 1154.1646 # Pa
    T_e = 1510.0513 # K
    rho_e = 0.0014269262 # kg/m^3
    m_dot_e = 4.061754 * n_engines #kg/s
    gamma = 1.2250769
    Ma = 4.3094718

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'bluemoon_H60_lambda1p5_m3p75e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 37500 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 0.694 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 60 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 12.1716 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 9.859006 # m/s (downward)
    a_thrust = 0.81 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  1714083 # Pa
    T_0 = 3376.0562 # K
    rho_0 = 0.86858998 # kg/m3

    # List Exit Conditions
    P_e = 1887.3576 # Pa
    T_e = 1490.7043 # K
    rho_e = 0.0023636825 # kg/m^3
    m_dot_e = 6.7544524 * n_engines #kg/s
    gamma = 1.2266141
    Ma = 4.3480774

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)



if spacecraft.lander_type == 'bluemoon_H45_lambda1p5_m1p5e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 15000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 0.694 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 45 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 10.5409 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 8.5381497 # m/s (downward)
    a_thrust = 0.81 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  682350.16 # Pa
    T_0 = 3256.3053 # K
    rho_0 = 0.35422271 # kg/m3

    # List Exit Conditions
    P_e = 781.94578 # Pa
    T_e = 1526.5296 # K
    rho_e = 0.00095631794 # kg/m^3
    m_dot_e = 2.7130224 * n_engines #kg/s
    gamma = 1.2235941
    Ma = 4.2774114

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'bluemoon_H90_lambda1p5_m1p5e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 15000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 0.694 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 90 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 14.9071 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 12.074767 # m/s (downward)
    a_thrust = 0.81 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  682350.16 # Pa
    T_0 = 3256.3053 # K
    rho_0 = 0.35422271 # kg/m3

    # List Exit Conditions
    P_e = 781.94578 # Pa
    T_e = 1526.5296 # K
    rho_e = 0.00095631794 # kg/m^3
    m_dot_e = 2.7130224 * n_engines #kg/s
    gamma = 1.2235941
    Ma = 4.2774114

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'bluemoon_H45_lambda1p5_m2p25e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 22500 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 0.694 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 45 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 10.5409 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 8.5381497 # m/s (downward)
    a_thrust = 0.81 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  1025641.4 # Pa
    T_0 = 3309.1795 # K
    rho_0 = 0.52670559 # kg/m3

    # List Exit Conditions
    P_e = 1154.1646 # Pa
    T_e = 1510.0513 # K
    rho_e = 0.0014269262 # kg/m^3
    m_dot_e = 4.061754 * n_engines #kg/s
    gamma = 1.2250769
    Ma = 4.3094718

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'bluemoon_H90_lambda1p5_m2p25e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 22500 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 0.694 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 90 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 14.9071 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 12.074767 # m/s (downward)
    a_thrust = 0.81 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  1025641.4 # Pa
    T_0 = 3309.1795 # K
    rho_0 = 0.52670559 # kg/m3

    # List Exit Conditions
    P_e = 1154.1646 # Pa
    T_e = 1510.0513 # K
    rho_e = 0.0014269262 # kg/m^3
    m_dot_e = 4.061754 * n_engines #kg/s
    gamma = 1.2250769
    Ma = 4.3094718

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'bluemoon_H45_lambda1p5_m3e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 30000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 0.694 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 45 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 10.5409 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 8.5381497 # m/s (downward)
    a_thrust = 0.81 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  1369629.5 # Pa
    T_0 = 3346.7844 # K
    rho_0 = 0.69806498 # kg/m3

    # List Exit Conditions
    P_e = 1522.2554 # Pa
    T_e = 1498.9504 # K
    rho_e = 0.0018958658 # kg/m^3
    m_dot_e = 5.4086981 * n_engines #kg/s
    gamma = 1.2259393
    Ma = 4.3314993

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'bluemoon_H90_lambda1p5_m3e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 30000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 0.694 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 90 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 14.9071 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 12.074767 # m/s (downward)
    a_thrust = 0.81 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  1369629.5 # Pa
    T_0 = 3346.7844 # K
    rho_0 = 0.69806498 # kg/m3

    # List Exit Conditions
    P_e = 1522.2554 # Pa
    T_e = 1498.9504 # K
    rho_e = 0.0018958658 # kg/m^3
    m_dot_e = 5.4086981 * n_engines #kg/s
    gamma = 1.2259393
    Ma = 4.3314993

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'bluemoon_H45_lambda1p5_m3p75e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 37500 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 0.694 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 45 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 10.5409 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 8.5381497 # m/s (downward)
    a_thrust = 0.81 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  1714083 # Pa
    T_0 = 3376.0562 # K
    rho_0 = 0.86858998 # kg/m3

    # List Exit Conditions
    P_e = 1887.3576 # Pa
    T_e = 1490.7043 # K
    rho_e = 0.0023636825 # kg/m^3
    m_dot_e = 6.7544524 * n_engines #kg/s
    gamma = 1.2266141
    Ma = 4.3480774

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'bluemoon_H90_lambda1p5_m3p75e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 37500 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 0.694 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 90 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 14.9071 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 12.074767 # m/s (downward)
    a_thrust = 0.81 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  1714083 # Pa
    T_0 = 3376.0562 # K
    rho_0 = 0.86858998 # kg/m3

    # List Exit Conditions
    P_e = 1887.3576 # Pa
    T_e = 1490.7043 # K
    rho_e = 0.0023636825 # kg/m^3
    m_dot_e = 6.7544524 * n_engines #kg/s
    gamma = 1.2266141
    Ma = 4.3480774

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'bluemoon_H45_lambda1p5_m4p5e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 45000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 0.694 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 45 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 10.5409 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 8.5381497 # m/s (downward)
    a_thrust = 0.81 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  2058928 # Pa
    T_0 = 3400.0125 # K
    rho_0 = 1.0384926 # kg/m3

    # List Exit Conditions
    P_e = 2250.2155 # Pa
    T_e = 1484.0982 # K
    rho_e = 0.0028305296 # kg/m^3
    m_dot_e = 8.0993142 * n_engines #kg/s
    gamma = 1.2271179
    Ma = 4.3612625

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'bluemoon_H90_lambda1p5_m4p5e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 45000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 0.694 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 90 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 14.9071 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 12.074767 # m/s (downward)
    a_thrust = 0.81 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  2058928 # Pa
    T_0 = 3400.0125 # K
    rho_0 = 1.0384926 # kg/m3

    # List Exit Conditions
    P_e = 2250.2155 # Pa
    T_e = 1484.0982 # K
    rho_e = 0.0028305296 # kg/m^3
    m_dot_e = 8.0993142 * n_engines #kg/s
    gamma = 1.2271179
    Ma = 4.3612625

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)



if spacecraft.lander_type == 'bluemoon_H45_lambda2_m1p5e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 15000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 0.694 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 45 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 7.45356 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 12.074767 # m/s (downward)
    a_thrust = 1.62 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  911141.54 # Pa
    T_0 = 3293.8264 # K
    rho_0 = 0.46937191 # kg/m3

    # List Exit Conditions
    P_e = 1030.7217 # Pa
    T_e = 1514.7321 # K
    rho_e = 0.0012702528 # kg/m^3
    m_dot_e = 3.612218 * n_engines #kg/s
    gamma = 1.2246446
    Ma = 4.3003136

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'bluemoon_H45_lambda2p5_m1p5e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 15000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 0.694 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 45 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 6.08581 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 14.788509 # m/s (downward)
    a_thrust = 2.43 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  1140249.8 # Pa
    T_0 = 3322.9482 # K
    rho_0 = 0.58393856 # kg/m3

    # List Exit Conditions
    P_e = 1277.2478 # Pa
    T_e = 1505.9718 # K
    rho_e = 0.0015834008 # kg/m^3
    m_dot_e = 4.5109027 * n_engines #kg/s
    gamma = 1.2254207
    Ma = 4.3175967

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'bluemoon_H45_lambda3_m1p5e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 15000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 0.694 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 45 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 5.27046 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 17.076299 # m/s (downward)
    a_thrust = 3.24 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  1369629.5 # Pa
    T_0 = 3346.7844 # K
    rho_0 = 0.69806498 # kg/m3

    # List Exit Conditions
    P_e = 1522.2554 # Pa
    T_e = 1498.9504 # K
    rho_e = 0.0018958658 # kg/m^3
    m_dot_e = 5.4086981 * n_engines #kg/s
    gamma = 1.2259393
    Ma = 4.3314993

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'bluemoon_H45_lambda4_m1p5e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 15000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 0.694 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 45 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 4.30331 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 20.91411 # m/s (downward)
    a_thrust = 4.86 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  1829007.5 # Pa
    T_0 = 3384.5595 # K
    rho_0 = 0.92528781 # kg/m3

    # List Exit Conditions
    P_e = 2008.5449 # Pa
    T_e = 1488.3198 # K
    rho_e = 0.0025193762 # kg/m^3
    m_dot_e = 7.2029886 * n_engines #kg/s
    gamma = 1.226758
    Ma = 4.3527313

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'bluemoon_H45_lambda5_m1p5e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 15000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 0.694 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 45 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 3.72678 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 24.149534 # m/s (downward)
    a_thrust = 6.48 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  2289051.9 # Pa
    T_0 = 3413.865 # K
    rho_0 = 1.1514354 # kg/m3

    # List Exit Conditions
    P_e = 2491.1483 # Pa
    T_e = 1480.3752 # K
    rho_e = 0.0031414201 # kg/m^3
    m_dot_e = 8.9954512 * n_engines #kg/s
    gamma = 1.2273891
    Ma = 4.3687606

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'bluemoon_H90_lambda2_m1p5e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 15000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 0.694 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 90 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 10.5409 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 17.076299 # m/s (downward)
    a_thrust = 1.62 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  911141.54 # Pa
    T_0 = 3293.8264 # K
    rho_0 = 0.46937191 # kg/m3

    # List Exit Conditions
    P_e = 1030.7217 # Pa
    T_e = 1514.7321 # K
    rho_e = 0.0012702528 # kg/m^3
    m_dot_e = 3.612218 * n_engines #kg/s
    gamma = 1.2246446
    Ma = 4.3003136

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'bluemoon_H90_lambda2p5_m1p5e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 15000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 0.694 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 90 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 8.60663 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 20.91411 # m/s (downward)
    a_thrust = 2.43 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  1140249.8 # Pa
    T_0 = 3322.9482 # K
    rho_0 = 0.58393856 # kg/m3

    # List Exit Conditions
    P_e = 1277.2478 # Pa
    T_e = 1505.9718 # K
    rho_e = 0.0015834008 # kg/m^3
    m_dot_e = 4.5109027 * n_engines #kg/s
    gamma = 1.2254207
    Ma = 4.3175967

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'bluemoon_H90_lambda3_m1p5e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 15000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 0.694 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 90 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 7.45356 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 24.149534 # m/s (downward)
    a_thrust = 3.24 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  1369629.5 # Pa
    T_0 = 3346.7844 # K
    rho_0 = 0.69806498 # kg/m3

    # List Exit Conditions
    P_e = 1522.2554 # Pa
    T_e = 1498.9504 # K
    rho_e = 0.0018958658 # kg/m^3
    m_dot_e = 5.4086981 * n_engines #kg/s
    gamma = 1.2259393
    Ma = 4.3314993

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'bluemoon_H90_lambda4_m1p5e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 15000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 0.694 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 90 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 6.08581 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 29.577018 # m/s (downward)
    a_thrust = 4.86 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  1829007.5 # Pa
    T_0 = 3384.5595 # K
    rho_0 = 0.92528781 # kg/m3

    # List Exit Conditions
    P_e = 2008.5449 # Pa
    T_e = 1488.3198 # K
    rho_e = 0.0025193762 # kg/m^3
    m_dot_e = 7.2029886 * n_engines #kg/s
    gamma = 1.226758
    Ma = 4.3527313

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'bluemoon_H90_lambda5_m1p5e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 15000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 0.694 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 90 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 5.27046 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 34.152599 # m/s (downward)
    a_thrust = 6.48 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  2289051.9 # Pa
    T_0 = 3413.865 # K
    rho_0 = 1.1514354 # kg/m3

    # List Exit Conditions
    P_e = 2491.1483 # Pa
    T_e = 1480.3752 # K
    rho_e = 0.0031414201 # kg/m^3
    m_dot_e = 8.9954512 * n_engines #kg/s
    gamma = 1.2273891
    Ma = 4.3687606

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)


# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------

if spacecraft.lander_type == 'starship_H30_lambda1p5_m5e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 50000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 5.89646 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 30 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 8.60663 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 6.97137 # m/s (downward)
    a_thrust = 0.81 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  366895.05 # Pa
    T_0 = 3057.554 # K
    rho_0 = 0.28819245 # kg/m3

    # List Exit Conditions
    P_e = 316.86731 # Pa
    T_e = 1827.0265 # K
    rho_e = 0.00053739418 # kg/m^3
    m_dot_e = 11.372391 * n_engines #kg/s
    gamma = 1.145966
    Ma = 4.1614505

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'starship_H30_lambda2_m5e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 50000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 5.89646 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 30 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 6.08581 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 9.859006 # m/s (downward)
    a_thrust = 1.62 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  489550.95 # Pa
    T_0 = 3077.9836 # K
    rho_0 = 0.38203972 # kg/m3

    # List Exit Conditions
    P_e = 418.13099 # Pa
    T_e = 1820.4668 # K
    rho_e = 0.00071343288 # kg/m^3
    m_dot_e = 15.167359 * n_engines #kg/s
    gamma = 1.148478
    Ma = 4.1714004

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'starship_H30_lambda4_m5e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 50000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 5.89646 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 30 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 3.51364 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 17.076299 # m/s (downward)
    a_thrust = 4.86 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  980174.57 # Pa
    T_0 = 3159.7019 # K
    rho_0 = 0.7574288 # kg/m3

    # List Exit Conditions
    P_e = 823.18573 # Pa
    T_e = 1794.2283 # K
    rho_e = 0.0014175877 # kg/m^3
    m_dot_e = 30.049521 * n_engines #kg/s
    gamma = 1.158526
    Ma = 4.2111998

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'starship_H30_lambda6_m5e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 50000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 5.89646 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 30 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 2.72166 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 22.045408 # m/s (downward)
    a_thrust = 8.1 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  1470798.2 # Pa
    T_0 = 3241.4201 # K
    rho_0 = 1.1328179 # kg/m3

    # List Exit Conditions
    P_e = 1228.2405 # Pa
    T_e = 1767.9897 # K
    rho_e = 0.0021217424 # kg/m^3
    m_dot_e = 44.46373 * n_engines #kg/s
    gamma = 1.1685739
    Ma = 4.2509991

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'starship_H30_lambda8_m5e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 50000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 5.89646 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 30 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 2.30022 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 26.084478 # m/s (downward)
    a_thrust = 11.34 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  1961421.8 # Pa
    T_0 = 3323.1384 # K
    rho_0 = 1.508207 # kg/m3

    # List Exit Conditions
    P_e = 1633.2952 # Pa
    T_e = 1741.7512 # K
    rho_e = 0.0028258972 # kg/m^3
    m_dot_e = 58.425179 * n_engines #kg/s
    gamma = 1.1786219
    Ma = 4.2907985

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'starship_H30_lambda10_m5e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 50000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 5.89646 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 30 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 2.0286 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 29.577018 # m/s (downward)
    a_thrust = 14.58 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  2452045.4 # Pa
    T_0 = 3404.8567 # K
    rho_0 = 1.883596 # kg/m3

    # List Exit Conditions
    P_e = 2038.3499 # Pa
    T_e = 1715.5126 # K
    rho_e = 0.003530052 # kg/m^3
    m_dot_e = 71.949251 * n_engines #kg/s
    gamma = 1.1886699
    Ma = 4.3305979

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'starship_H60_lambda1p5_m5e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 50000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 5.89646 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 60 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 12.1716 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 9.859006 # m/s (downward)
    a_thrust = 0.81 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  366895.05 # Pa
    T_0 = 3057.554 # K
    rho_0 = 0.28819245 # kg/m3

    # List Exit Conditions
    P_e = 316.86731 # Pa
    T_e = 1827.0265 # K
    rho_e = 0.00053739418 # kg/m^3
    m_dot_e = 11.372391 * n_engines #kg/s
    gamma = 1.145966
    Ma = 4.1614505

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'starship_H60_lambda2_m5e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 50000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 5.89646 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 60 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 8.60663 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 13.94274 # m/s (downward)
    a_thrust = 1.62 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  489550.95 # Pa
    T_0 = 3077.9836 # K
    rho_0 = 0.38203972 # kg/m3

    # List Exit Conditions
    P_e = 418.13099 # Pa
    T_e = 1820.4668 # K
    rho_e = 0.00071343288 # kg/m^3
    m_dot_e = 15.167359 * n_engines #kg/s
    gamma = 1.148478
    Ma = 4.1714004

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'starship_H60_lambda4_m5e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 50000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 5.89646 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 60 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 4.96904 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 24.149534 # m/s (downward)
    a_thrust = 4.86 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  980174.57 # Pa
    T_0 = 3159.7019 # K
    rho_0 = 0.7574288 # kg/m3

    # List Exit Conditions
    P_e = 823.18573 # Pa
    T_e = 1794.2283 # K
    rho_e = 0.0014175877 # kg/m^3
    m_dot_e = 30.049521 * n_engines #kg/s
    gamma = 1.158526
    Ma = 4.2111998

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'starship_H60_lambda6_m5e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 50000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 5.89646 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 60 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 3.849 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 31.176915 # m/s (downward)
    a_thrust = 8.1 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  1470798.2 # Pa
    T_0 = 3241.4201 # K
    rho_0 = 1.1328179 # kg/m3

    # List Exit Conditions
    P_e = 1228.2405 # Pa
    T_e = 1767.9897 # K
    rho_e = 0.0021217424 # kg/m^3
    m_dot_e = 44.46373 * n_engines #kg/s
    gamma = 1.1685739
    Ma = 4.2509991

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'starship_H60_lambda8_m5e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 50000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 5.89646 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 60 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 3.253 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 36.889023 # m/s (downward)
    a_thrust = 11.34 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  1961421.8 # Pa
    T_0 = 3323.1384 # K
    rho_0 = 1.508207 # kg/m3

    # List Exit Conditions
    P_e = 1633.2952 # Pa
    T_e = 1741.7512 # K
    rho_e = 0.0028258972 # kg/m^3
    m_dot_e = 58.425179 * n_engines #kg/s
    gamma = 1.1786219
    Ma = 4.2907985

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'starship_H60_lambda10_m5e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 50000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 5.89646 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 60 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 2.86888 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 41.82822 # m/s (downward)
    a_thrust = 14.58 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  2452045.4 # Pa
    T_0 = 3404.8567 # K
    rho_0 = 1.883596 # kg/m3

    # List Exit Conditions
    P_e = 2038.3499 # Pa
    T_e = 1715.5126 # K
    rho_e = 0.003530052 # kg/m^3
    m_dot_e = 71.949251 * n_engines #kg/s
    gamma = 1.1886699
    Ma = 4.3305979

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'starship_H30_lambda1p5_m1e05':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 100000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 5.89646 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 30 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 8.60663 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 6.97137 # m/s (downward)
    a_thrust = 0.81 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  734862.76 # Pa
    T_0 = 3118.8427 # K
    rho_0 = 0.56973426 # kg/m3

    # List Exit Conditions
    P_e = 620.65836 # Pa
    T_e = 1807.3475 # K
    rho_e = 0.0010655103 # kg/m^3
    m_dot_e = 22.667799 * n_engines #kg/s
    gamma = 1.153502
    Ma = 4.1913001

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'starship_H30_lambda2_m1e05':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 100000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 5.89646 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 30 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 6.08581 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 9.859006 # m/s (downward)
    a_thrust = 1.62 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  980174.57 # Pa
    T_0 = 3159.7019 # K
    rho_0 = 0.7574288 # kg/m3

    # List Exit Conditions
    P_e = 823.18573 # Pa
    T_e = 1794.2283 # K
    rho_e = 0.0014175877 # kg/m^3
    m_dot_e = 30.049521 * n_engines #kg/s
    gamma = 1.158526
    Ma = 4.2111998

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'starship_H30_lambda4_m1e05':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 100000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 5.89646 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 30 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 3.51364 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 17.076299 # m/s (downward)
    a_thrust = 4.86 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  1961421.8 # Pa
    T_0 = 3323.1384 # K
    rho_0 = 1.508207 # kg/m3

    # List Exit Conditions
    P_e = 1633.2952 # Pa
    T_e = 1741.7512 # K
    rho_e = 0.0028258972 # kg/m^3
    m_dot_e = 58.425179 * n_engines #kg/s
    gamma = 1.1786219
    Ma = 4.2907985

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'starship_H30_lambda6_m1e05':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 100000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 5.89646 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 30 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 2.72166 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 22.045408 # m/s (downward)
    a_thrust = 8.1 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  2944532.7 # Pa
    T_0 = 3442.0548 # K
    rho_0 = 2.2521125 # kg/m3

    # List Exit Conditions
    P_e = 2429.5801 # Pa
    T_e = 1703.5621 # K
    rho_e = 0.0042257519 # kg/m^3
    m_dot_e = 86.015667 * n_engines #kg/s
    gamma = 1.1920032
    Ma = 4.3501659

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'starship_H30_lambda8_m1e05':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 100000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 5.89646 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 30 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 2.30022 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 26.084478 # m/s (downward)
    a_thrust = 11.34 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  3931116.7 # Pa
    T_0 = 3478.0059 # K
    rho_0 = 2.9832108 # kg/m3

    # List Exit Conditions
    P_e = 3200.1022 # Pa
    T_e = 1691.9993 # K
    rho_e = 0.0056098503 # kg/m^3
    m_dot_e = 114.61105 * n_engines #kg/s
    gamma = 1.1928714
    Ma = 4.3718313

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'starship_H30_lambda10_m1e05':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 100000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 5.89646 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 30 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 2.0286 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 29.577018 # m/s (downward)
    a_thrust = 14.58 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  4917700.8 # Pa
    T_0 = 3513.957 # K
    rho_0 = 3.714309 # kg/m3

    # List Exit Conditions
    P_e = 3970.6243 # Pa
    T_e = 1680.4365 # K
    rho_e = 0.0069939488 # kg/m^3
    m_dot_e = 142.84175 * n_engines #kg/s
    gamma = 1.1937396
    Ma = 4.3934967

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'starship_H60_lambda1p5_m1e05':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 100000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 5.89646 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 60 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 12.1716 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 9.859006 # m/s (downward)
    a_thrust = 0.81 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  734862.76 # Pa
    T_0 = 3118.8427 # K
    rho_0 = 0.56973426 # kg/m3

    # List Exit Conditions
    P_e = 620.65836 # Pa
    T_e = 1807.3475 # K
    rho_e = 0.0010655103 # kg/m^3
    m_dot_e = 22.667799 * n_engines #kg/s
    gamma = 1.153502
    Ma = 4.1913001

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'starship_H60_lambda2_m1e05':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 100000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 5.89646 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 60 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 8.60663 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 13.94274 # m/s (downward)
    a_thrust = 1.62 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  980174.57 # Pa
    T_0 = 3159.7019 # K
    rho_0 = 0.7574288 # kg/m3

    # List Exit Conditions
    P_e = 823.18573 # Pa
    T_e = 1794.2283 # K
    rho_e = 0.0014175877 # kg/m^3
    m_dot_e = 30.049521 * n_engines #kg/s
    gamma = 1.158526
    Ma = 4.2111998

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'starship_H60_lambda4_m1e05':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 100000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 5.89646 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 60 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 4.96904 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 24.149534 # m/s (downward)
    a_thrust = 4.86 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  1961421.8 # Pa
    T_0 = 3323.1384 # K
    rho_0 = 1.508207 # kg/m3

    # List Exit Conditions
    P_e = 1633.2952 # Pa
    T_e = 1741.7512 # K
    rho_e = 0.0028258972 # kg/m^3
    m_dot_e = 58.425179 * n_engines #kg/s
    gamma = 1.1786219
    Ma = 4.2907985

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'starship_H60_lambda6_m1e05':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 100000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 5.89646 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 60 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 3.849 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 31.176915 # m/s (downward)
    a_thrust = 8.1 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  2944532.7 # Pa
    T_0 = 3442.0548 # K
    rho_0 = 2.2521125 # kg/m3

    # List Exit Conditions
    P_e = 2429.5801 # Pa
    T_e = 1703.5621 # K
    rho_e = 0.0042257519 # kg/m^3
    m_dot_e = 86.015667 * n_engines #kg/s
    gamma = 1.1920032
    Ma = 4.3501659

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'starship_H60_lambda8_m1e05':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 100000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 5.89646 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 60 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 3.253 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 36.889023 # m/s (downward)
    a_thrust = 11.34 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  3931116.7 # Pa
    T_0 = 3478.0059 # K
    rho_0 = 2.9832108 # kg/m3

    # List Exit Conditions
    P_e = 3200.1022 # Pa
    T_e = 1691.9993 # K
    rho_e = 0.0056098503 # kg/m^3
    m_dot_e = 114.61105 * n_engines #kg/s
    gamma = 1.1928714
    Ma = 4.3718313

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'starship_H60_lambda10_m1e05':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 100000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 5.89646 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 60 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 1 # m
    lander_total_sim_time = 2.86888 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 41.82822 # m/s (downward)
    a_thrust = 14.58 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  4917700.8 # Pa
    T_0 = 3513.957 # K
    rho_0 = 3.714309 # kg/m3

    # List Exit Conditions
    P_e = 3970.6243 # Pa
    T_e = 1680.4365 # K
    rho_e = 0.0069939488 # kg/m^3
    m_dot_e = 142.84175 * n_engines #kg/s
    gamma = 1.1937396
    Ma = 4.3934967

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

# ------------------------------------------------------------------------------------------------
# 11/17/2025
# ------------------------------------------------------------------------------------------------

if spacecraft.lander_type == 'starshipUpper_H30_lambda1p5_m5e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 50000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 5.89646 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 30 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 33 # m
    lander_total_sim_time = 8.60663 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 6.97137 # m/s (downward)
    a_thrust = 0.81 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  366895.05 # Pa
    T_0 = 3057.554 # K
    rho_0 = 0.28819245 # kg/m3

    # List Exit Conditions
    P_e = 316.86731 # Pa
    T_e = 1827.0265 # K
    rho_e = 0.00053739418 # kg/m^3
    m_dot_e = 11.372391 * n_engines #kg/s
    gamma = 1.145966
    Ma = 4.1614505

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'starshipUpper_H30_lambda2_m5e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 50000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 5.89646 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 30 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 33 # m
    lander_total_sim_time = 6.08581 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 9.859006 # m/s (downward)
    a_thrust = 1.62 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  489550.95 # Pa
    T_0 = 3077.9836 # K
    rho_0 = 0.38203972 # kg/m3

    # List Exit Conditions
    P_e = 418.13099 # Pa
    T_e = 1820.4668 # K
    rho_e = 0.00071343288 # kg/m^3
    m_dot_e = 15.167359 * n_engines #kg/s
    gamma = 1.148478
    Ma = 4.1714004

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'starshipUpper_H30_lambda4_m5e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 50000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 5.89646 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 30 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 33 # m
    lander_total_sim_time = 3.51364 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 17.076299 # m/s (downward)
    a_thrust = 4.86 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  980174.57 # Pa
    T_0 = 3159.7019 # K
    rho_0 = 0.7574288 # kg/m3

    # List Exit Conditions
    P_e = 823.18573 # Pa
    T_e = 1794.2283 # K
    rho_e = 0.0014175877 # kg/m^3
    m_dot_e = 30.049521 * n_engines #kg/s
    gamma = 1.158526
    Ma = 4.2111998

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'starshipUpper_H30_lambda6_m5e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 50000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 5.89646 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 30 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 33 # m
    lander_total_sim_time = 2.72166 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 22.045408 # m/s (downward)
    a_thrust = 8.1 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  1470798.2 # Pa
    T_0 = 3241.4201 # K
    rho_0 = 1.1328179 # kg/m3

    # List Exit Conditions
    P_e = 1228.2405 # Pa
    T_e = 1767.9897 # K
    rho_e = 0.0021217424 # kg/m^3
    m_dot_e = 44.46373 * n_engines #kg/s
    gamma = 1.1685739
    Ma = 4.2509991

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'starshipUpper_H30_lambda8_m5e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 50000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 5.89646 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 30 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 33 # m
    lander_total_sim_time = 2.30022 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 26.084478 # m/s (downward)
    a_thrust = 11.34 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  1961421.8 # Pa
    T_0 = 3323.1384 # K
    rho_0 = 1.508207 # kg/m3

    # List Exit Conditions
    P_e = 1633.2952 # Pa
    T_e = 1741.7512 # K
    rho_e = 0.0028258972 # kg/m^3
    m_dot_e = 58.425179 * n_engines #kg/s
    gamma = 1.1786219
    Ma = 4.2907985

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'starshipUpper_H30_lambda10_m5e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 50000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 5.89646 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 30 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 33 # m
    lander_total_sim_time = 2.0286 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 29.577018 # m/s (downward)
    a_thrust = 14.58 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  2452045.4 # Pa
    T_0 = 3404.8567 # K
    rho_0 = 1.883596 # kg/m3

    # List Exit Conditions
    P_e = 2038.3499 # Pa
    T_e = 1715.5126 # K
    rho_e = 0.003530052 # kg/m^3
    m_dot_e = 71.949251 * n_engines #kg/s
    gamma = 1.1886699
    Ma = 4.3305979

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'starshipUpper_H60_lambda1p5_m5e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 50000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 5.89646 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 60 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 33 # m
    lander_total_sim_time = 12.1716 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 9.859006 # m/s (downward)
    a_thrust = 0.81 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  366895.05 # Pa
    T_0 = 3057.554 # K
    rho_0 = 0.28819245 # kg/m3

    # List Exit Conditions
    P_e = 316.86731 # Pa
    T_e = 1827.0265 # K
    rho_e = 0.00053739418 # kg/m^3
    m_dot_e = 11.372391 * n_engines #kg/s
    gamma = 1.145966
    Ma = 4.1614505

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'starshipUpper_H60_lambda2_m5e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 50000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 5.89646 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 60 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 33 # m
    lander_total_sim_time = 8.60663 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 13.94274 # m/s (downward)
    a_thrust = 1.62 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  489550.95 # Pa
    T_0 = 3077.9836 # K
    rho_0 = 0.38203972 # kg/m3

    # List Exit Conditions
    P_e = 418.13099 # Pa
    T_e = 1820.4668 # K
    rho_e = 0.00071343288 # kg/m^3
    m_dot_e = 15.167359 * n_engines #kg/s
    gamma = 1.148478
    Ma = 4.1714004

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'starshipUpper_H60_lambda4_m5e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 50000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 5.89646 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 60 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 33 # m
    lander_total_sim_time = 4.96904 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 24.149534 # m/s (downward)
    a_thrust = 4.86 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  980174.57 # Pa
    T_0 = 3159.7019 # K
    rho_0 = 0.7574288 # kg/m3

    # List Exit Conditions
    P_e = 823.18573 # Pa
    T_e = 1794.2283 # K
    rho_e = 0.0014175877 # kg/m^3
    m_dot_e = 30.049521 * n_engines #kg/s
    gamma = 1.158526
    Ma = 4.2111998

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'starshipUpper_H60_lambda6_m5e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 50000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 5.89646 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 60 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 33 # m
    lander_total_sim_time = 3.849 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 31.176915 # m/s (downward)
    a_thrust = 8.1 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  1470798.2 # Pa
    T_0 = 3241.4201 # K
    rho_0 = 1.1328179 # kg/m3

    # List Exit Conditions
    P_e = 1228.2405 # Pa
    T_e = 1767.9897 # K
    rho_e = 0.0021217424 # kg/m^3
    m_dot_e = 44.46373 * n_engines #kg/s
    gamma = 1.1685739
    Ma = 4.2509991

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'starshipUpper_H60_lambda8_m5e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 50000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 5.89646 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 60 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 33 # m
    lander_total_sim_time = 3.253 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 36.889023 # m/s (downward)
    a_thrust = 11.34 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  1961421.8 # Pa
    T_0 = 3323.1384 # K
    rho_0 = 1.508207 # kg/m3

    # List Exit Conditions
    P_e = 1633.2952 # Pa
    T_e = 1741.7512 # K
    rho_e = 0.0028258972 # kg/m^3
    m_dot_e = 58.425179 * n_engines #kg/s
    gamma = 1.1786219
    Ma = 4.2907985

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'starshipUpper_H60_lambda10_m5e04':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 50000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 5.89646 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 60 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 33 # m
    lander_total_sim_time = 2.86888 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 41.82822 # m/s (downward)
    a_thrust = 14.58 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  2452045.4 # Pa
    T_0 = 3404.8567 # K
    rho_0 = 1.883596 # kg/m3

    # List Exit Conditions
    P_e = 2038.3499 # Pa
    T_e = 1715.5126 # K
    rho_e = 0.003530052 # kg/m^3
    m_dot_e = 71.949251 * n_engines #kg/s
    gamma = 1.1886699
    Ma = 4.3305979

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'starshipUpper_H30_lambda1p5_m1e05':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 100000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 5.89646 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 30 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 33 # m
    lander_total_sim_time = 8.60663 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 6.97137 # m/s (downward)
    a_thrust = 0.81 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  734862.76 # Pa
    T_0 = 3118.8427 # K
    rho_0 = 0.56973426 # kg/m3

    # List Exit Conditions
    P_e = 620.65836 # Pa
    T_e = 1807.3475 # K
    rho_e = 0.0010655103 # kg/m^3
    m_dot_e = 22.667799 * n_engines #kg/s
    gamma = 1.153502
    Ma = 4.1913001

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'starshipUpper_H30_lambda2_m1e05':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 100000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 5.89646 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 30 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 33 # m
    lander_total_sim_time = 6.08581 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 9.859006 # m/s (downward)
    a_thrust = 1.62 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  980174.57 # Pa
    T_0 = 3159.7019 # K
    rho_0 = 0.7574288 # kg/m3

    # List Exit Conditions
    P_e = 823.18573 # Pa
    T_e = 1794.2283 # K
    rho_e = 0.0014175877 # kg/m^3
    m_dot_e = 30.049521 * n_engines #kg/s
    gamma = 1.158526
    Ma = 4.2111998

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'starshipUpper_H30_lambda4_m1e05':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 100000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 5.89646 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 30 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 33 # m
    lander_total_sim_time = 3.51364 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 17.076299 # m/s (downward)
    a_thrust = 4.86 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  1961421.8 # Pa
    T_0 = 3323.1384 # K
    rho_0 = 1.508207 # kg/m3

    # List Exit Conditions
    P_e = 1633.2952 # Pa
    T_e = 1741.7512 # K
    rho_e = 0.0028258972 # kg/m^3
    m_dot_e = 58.425179 * n_engines #kg/s
    gamma = 1.1786219
    Ma = 4.2907985

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'starshipUpper_H30_lambda6_m1e05':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 100000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 5.89646 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 30 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 33 # m
    lander_total_sim_time = 2.72166 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 22.045408 # m/s (downward)
    a_thrust = 8.1 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  2944532.7 # Pa
    T_0 = 3442.0548 # K
    rho_0 = 2.2521125 # kg/m3

    # List Exit Conditions
    P_e = 2429.5801 # Pa
    T_e = 1703.5621 # K
    rho_e = 0.0042257519 # kg/m^3
    m_dot_e = 86.015667 * n_engines #kg/s
    gamma = 1.1920032
    Ma = 4.3501659

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'starshipUpper_H30_lambda8_m1e05':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 100000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 5.89646 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 30 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 33 # m
    lander_total_sim_time = 2.30022 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 26.084478 # m/s (downward)
    a_thrust = 11.34 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  3931116.7 # Pa
    T_0 = 3478.0059 # K
    rho_0 = 2.9832108 # kg/m3

    # List Exit Conditions
    P_e = 3200.1022 # Pa
    T_e = 1691.9993 # K
    rho_e = 0.0056098503 # kg/m^3
    m_dot_e = 114.61105 * n_engines #kg/s
    gamma = 1.1928714
    Ma = 4.3718313

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'starshipUpper_H30_lambda10_m1e05':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 100000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 5.89646 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 30 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 33 # m
    lander_total_sim_time = 2.0286 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 29.577018 # m/s (downward)
    a_thrust = 14.58 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  4917700.8 # Pa
    T_0 = 3513.957 # K
    rho_0 = 3.714309 # kg/m3

    # List Exit Conditions
    P_e = 3970.6243 # Pa
    T_e = 1680.4365 # K
    rho_e = 0.0069939488 # kg/m^3
    m_dot_e = 142.84175 * n_engines #kg/s
    gamma = 1.1937396
    Ma = 4.3934967

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'starshipUpper_H60_lambda1p5_m1e05':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 100000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 5.89646 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 60 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 33 # m
    lander_total_sim_time = 12.1716 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 9.859006 # m/s (downward)
    a_thrust = 0.81 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  734862.76 # Pa
    T_0 = 3118.8427 # K
    rho_0 = 0.56973426 # kg/m3

    # List Exit Conditions
    P_e = 620.65836 # Pa
    T_e = 1807.3475 # K
    rho_e = 0.0010655103 # kg/m^3
    m_dot_e = 22.667799 * n_engines #kg/s
    gamma = 1.153502
    Ma = 4.1913001

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'starshipUpper_H60_lambda2_m1e05':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 100000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 5.89646 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 60 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 33 # m
    lander_total_sim_time = 8.60663 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 13.94274 # m/s (downward)
    a_thrust = 1.62 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  980174.57 # Pa
    T_0 = 3159.7019 # K
    rho_0 = 0.7574288 # kg/m3

    # List Exit Conditions
    P_e = 823.18573 # Pa
    T_e = 1794.2283 # K
    rho_e = 0.0014175877 # kg/m^3
    m_dot_e = 30.049521 * n_engines #kg/s
    gamma = 1.158526
    Ma = 4.2111998

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'starshipUpper_H60_lambda4_m1e05':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 100000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 5.89646 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 60 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 33 # m
    lander_total_sim_time = 4.96904 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 24.149534 # m/s (downward)
    a_thrust = 4.86 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  1961421.8 # Pa
    T_0 = 3323.1384 # K
    rho_0 = 1.508207 # kg/m3

    # List Exit Conditions
    P_e = 1633.2952 # Pa
    T_e = 1741.7512 # K
    rho_e = 0.0028258972 # kg/m^3
    m_dot_e = 58.425179 * n_engines #kg/s
    gamma = 1.1786219
    Ma = 4.2907985

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'starshipUpper_H60_lambda6_m1e05':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 100000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 5.89646 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 60 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 33 # m
    lander_total_sim_time = 3.849 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 31.176915 # m/s (downward)
    a_thrust = 8.1 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  2944532.7 # Pa
    T_0 = 3442.0548 # K
    rho_0 = 2.2521125 # kg/m3

    # List Exit Conditions
    P_e = 2429.5801 # Pa
    T_e = 1703.5621 # K
    rho_e = 0.0042257519 # kg/m^3
    m_dot_e = 86.015667 * n_engines #kg/s
    gamma = 1.1920032
    Ma = 4.3501659

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'starshipUpper_H60_lambda8_m1e05':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 100000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 5.89646 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 60 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 33 # m
    lander_total_sim_time = 3.253 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 36.889023 # m/s (downward)
    a_thrust = 11.34 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  3931116.7 # Pa
    T_0 = 3478.0059 # K
    rho_0 = 2.9832108 # kg/m3

    # List Exit Conditions
    P_e = 3200.1022 # Pa
    T_e = 1691.9993 # K
    rho_e = 0.0056098503 # kg/m^3
    m_dot_e = 114.61105 * n_engines #kg/s
    gamma = 1.1928714
    Ma = 4.3718313

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

if spacecraft.lander_type == 'starshipUpper_H60_lambda10_m1e05':
    # Blue Origin uses three RE-7 Engines for its BlueMoon MK2 Lander
    n_engines = 3
    m_lander = 100000 # kg (Total Mass Lander)

    # List Thruster Exit Geometry
    A_nozzle = 5.89646 * n_engines  #m^2
    A_throat = A_nozzle/76 # m^2
    D_nozzle = 2 * np.sqrt(A_nozzle/np.pi)  #m^2
    r_nozzle = D_nozzle/2

    # Define Starting and Final Landing Positions, and Corresponding Nozzle Exit Pos
    lander_base_starting_altitude = 60 # m
    lander_base_final_altitude = 0 # m
    lander_base_to_nozzle_disp = 33 # m
    lander_total_sim_time = 2.86888 #
    h_nozzle_init = lander_base_starting_altitude + lander_base_to_nozzle_disp
    min_altitude_lander = lander_base_final_altitude + lander_base_to_nozzle_disp

    # Descent Properties (Downward Velocity and Upward Accelleration)
    v_init = 41.82822 # m/s (downward)
    a_thrust = 14.58 # m/s^2 (upward)

    # List Chamber Conditions
    P_0 =  4917700.8 # Pa
    T_0 = 3513.957 # K
    rho_0 = 3.714309 # kg/m3

    # List Exit Conditions
    P_e = 3970.6243 # Pa
    T_e = 1680.4365 # K
    rho_e = 0.0069939488 # kg/m^3
    m_dot_e = 142.84175 * n_engines #kg/s
    gamma = 1.1937396
    Ma = 4.3934967

    # Other Easily Computed Exit Variables
    R_gas = P_e / (rho_e * T_e)
    v_e = Ma * np.sqrt(gamma * R_gas * T_e)
    k_bar = gamma * (gamma - 1) * Ma**2
    timesteps.n_sub_timesteps = int(lander_total_sim_time/timesteps.delta_t)

import numpy as np
import var_range_of_interest as bounds

p_gas = None
T_gas = None
rho_gas = None
v_gas = None
v_T = None

rho_gas_arr = np.zeros(bounds.n_points_centerline)
v_gas_arr = np.zeros(bounds.n_points_centerline)
T_gas_arr = np.zeros(bounds.n_points_centerline)
v_T_arr = np.zeros(bounds.n_points_centerline)
p_gas_arr = np.zeros(bounds.n_points_centerline)    
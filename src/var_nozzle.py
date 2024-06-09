import numpy as np
from fun_nozzle import compute_gas_Mw
import var_constants as cs
from fun_compressible import *

# ---------------------------------------------------------
# Change only these variables (Input Conditions)
# ---------------------------------------------------------
gamma = 1.2
#P_0 = 19046 # Pascals (40 Ton)
P_0 = 69259.33 # Pascals (40 Ton)
T_0 = 4875 # Kelvin
Ma = 2.5
D_nozzle = 1.6 # meters
# ---------------------------------------------------------

M_w = compute_gas_Mw()

# determine gas/nozzle variables at the exit
R_gas = cs.R_universal/M_w # J/(kgK)

 # calculate stagnation density
rho_0 = P_0/(R_gas * T_0) # kg/m^3

 # calculate nozzle properties
A_nozzle = (np.pi/4) * D_nozzle**2 #m^2
A_throat = A_nozzle/A_over_Astar(Ma, gamma) # m^2

k_bar = gamma * (gamma - 1) * Ma**2


# initialize exit conditions
m_dot_e = None
P_e = None 
T_e = None 
rho_e = None
v_e = None
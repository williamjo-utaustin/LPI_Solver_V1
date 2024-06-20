import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate 
import sys

sys.path.append('/Users/williamjo/Documents/LPI/Codes/plume_regolith_solver_v1/src')

# import variables
import var_constants as cs
import var_nozzle as nozzle
import var_soil as soil
import var_impinged_gas as imp
import var_range_of_interest as bounds

# import functions
from fun_compressible import *
from fun_soil import 
from fun_nozzle import *
from fun_impinged_gas import *
from subroutines import *

# run main code
compute_nozzle_exhaust()
loop()

#plot_initial_disturbed_altitude()
#d_particle = 1
#
#def dudt(t,u_p):
#    C_d = 0.40
#    A_p = 4 * np.pi * d_particle**2
#    m_p = soil.rho_p * (4/3) * np.pi * (d_particle/2)**3
#    return (C_d * A_p * (rho_gas/2) * (v_gas - u_p)**2)/m_p
#
#t_span = np.linspace(0,1,1000)
#yinit = [0]
#sol = integrate.solve_ivp(lambda t, u_p: dudt(t,u_p),  [t_span[0], t_span[-1]], yinit, method='RK45')
#plt.semilogy(sol.t,sol.y[0])
#plt.show()
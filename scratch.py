import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate 
import sys
import time

sys.path.append('/Users/williamjo/Documents/LPI/Codes/plume_regolith_solver_v1/src')

# import variables
import var_constants as cs
import var_nozzle as nozzle
import var_soil as soil
import var_impinged_gas as imp
import var_range_of_interest as bounds

# import functions
from fun_compressible import *
from fun_soil import *
from fun_nozzle import *
from fun_impinged_gas import *
from subroutines import *


#ratio = (soil.rho_0 * cs.g * soil.D_bracket + soil.alpha_0)/ ((1100.700*1.62* soil.D_bracket + 0.2930747))
#print(soil.rho_0 * cs.g * soil.D_bracket + soil.alpha_0)
#print((1100.700*1.62* soil.D_bracket + 0.2930747))
#print(ratio)
#
#
param = np.asarray([multiplier * magnitude for magnitude in [1E-7, 1E-6, 1E-5, 1E-4, 1E-3, 1E-2] for multiplier in [1, 2, 3, 4, 5, 6, 7, 8, 9]])
#print(np.asarray(param))
#
g = 1.62      # Gravity m/s^2
A_const= 1.71575E-7
beta = 0.78
rho_p = 1100.35185

rho_gas = 6.929E-11
p_gas = 0.001376
T_gas = 4762
u_gas = 3665

k_b = 1.380649E-23
R = 406

mean_free_path = np.sqrt(np.pi/(2 * p_gas * rho_gas)) * A_const * T_gas**beta

def cd_sphere(Re, Kn):
    
    if Re <= 0.0:
        CD = 0.0
    elif Re > 8.0e6:
        CD = 0.2
    elif Re > 0.0 and Re <= 2:
        CD = (24.0/Re)
    elif Re > 2 and Re <= 100.0:
        p = np.array([4.22, -14.05, 34.87, 0.658])
        CD = np.polyval(p, 1.0/Re) 
    elif Re > 100.0 and Re <= 1.0e4:
        p = np.array([-30.41, 43.72, -17.08, 2.41])
        CD = np.polyval(p, 1.0/np.log10(Re))
    elif Re > 1.0e4 and Re <= 3.35e5:
        p = np.array([-0.1584, 2.031, -8.472, 11.932])
        CD = np.polyval(p, np.log10(Re))
    elif Re > 3.35e5 and Re <= 5.0e5:
        x1 = np.log10(Re/4.5e5)
        CD = 91.08*x1**4 + 0.0764
    else:
        p = np.array([-0.06338, 1.1905, -7.332, 14.93])
        CD = np.polyval(p, np.log10(Re))

    # add in correction factor to account for rarefied effects
    S_correction = 1 + Kn * (2.514 + 0.8 * np.exp(-0.55/Kn))
    CD = CD/S_correction
    return CD

def dudt(u, t, d_p):
    Kn_p = mean_free_path/d_p
    Re = (rho_gas * np.abs(u_gas - u) * d_p)/(A_const* T_gas ** beta)
    C_D = cd_sphere(Re, Kn_p)
    F_D = (3.0*rho_gas*C_D)*(u_gas - u)**2/(4.0*rho_p*d_p)
    return F_D + ((1 - (rho_gas/rho_p))*g)

def rk_4(func, u_n, t, dt, d_p):
    dt2 = dt/2.0
    k1 = func(u_n, t, d_p)
    k2 = func(u_n + k1*dt2, t + dt2, d_p)
    k3 = func(u_n + k2*dt2, t + dt2, d_p)
    k4 = func(u_n + k3*dt, t + dt, d_p)
    u_n_p1 = u_n + dt/6.0*(k1 + 2.0*k2 + 2.0*k3 + k4)
    return u_n_p1

dt = 0.001
n_timesteps = 100000
baseline_height = 0.3

u_ej = np.zeros_like(param)

for j in range(0,np.size(param)):
    
    d_p = param[j]
    t = 0
    u_p = 0

    x_n = 0
    for i in range(0,n_timesteps):
        t = t + dt 
        u_p_save = u_p
        u_p = rk_4(dudt, u_p, t, dt, d_p)

        x_n = x_n + u_p * dt
        print(i, t, x_n, u_p)
        if(x_n * np.tan(1.5*np.pi/180) > baseline_height or (u_p-u_p_save)/u_p < 0.00001):
            u_ej[j] = u_p
            break

plt.loglog(param, u_ej)
#plt.ylim(1,1E2)
plt.show()
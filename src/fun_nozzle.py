import numpy as np

import sys
sys.path.append('/Users/williamjo/Documents/LPI/Codes/plume_regolith_solver_v1/src')

import var_constants as cs
import var_nozzle as nozzle
import var_timestepping as timestep
from fun_compressible import *


def compute_gas_Mw():
    # insert dataset for gas molar mass
    gas_data = np.genfromtxt('/Users/williamjo/Documents/LPI/Codes/plume_regolith_solver_v1/data/FontesEtAl_GasData.csv',delimiter = ',')

    # compute molar mass from gas data
    M_w = 0
    for i in range(0,np.size(gas_data[:,0])):
        M_w = M_w + gas_data[i,0] * gas_data[i,1]

    return M_w/1000 # convert to kg/mol

def compute_nozzle_exhaust():
    nozzle.m_dot_e = mDot_over_A(nozzle.Ma, nozzle.gamma, nozzle.P_0, nozzle.R_gas, nozzle.T_0) * nozzle.A_nozzle
    nozzle.P_e = nozzle.P_0/(p0_over_p(nozzle.Ma, nozzle.gamma))
    nozzle.T_e = nozzle.T_0/(T0_over_T(nozzle.Ma, nozzle.gamma))
    nozzle.rho_e = (nozzle.P_0/(nozzle.R_gas*nozzle.T_0)) / (rho0_over_rho(nozzle.Ma, nozzle.gamma))
    nozzle.v_e = nozzle.Ma * np.sqrt(nozzle.gamma * nozzle.R_gas * nozzle.T_e)
    return None


def update_nozzle_height(h_nozzle):
    # update the new nozzle velocity
    h_nozzle = h_nozzle - (nozzle.v_descent * timestep.delta_t)
    return h_nozzle
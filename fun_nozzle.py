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

# WJ new notes
# 10/11/25
# set an upward accelleration to counteract the negative velocity
# need to make sure that the height will be 0 and the will be 0, even with displacement

def init_descent_velocity():
    nozzle.v_descent = nozzle.v_init

# we update the nozzle height to reflect the props at the next timestep (we put this routine in the current timestep)
def update_nozzle_height(h_nozzle, timestep_count):

    nozzle.v_descent = nozzle.v_init - nozzle.a_thrust * (timestep.delta_t * (timestep_count + 1))
    h_nozzle = nozzle.h_nozzle_init - nozzle.v_init * (timestep.delta_t * (timestep_count+1)) + 0.5 * nozzle.a_thrust * (timestep.delta_t * (timestep_count+1))**2
    
    return h_nozzle

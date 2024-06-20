import numpy as np
import matplotlib.pyplot as plt
import sys 

sys.path.append('/Users/williamjo/Documents/LPI/Codes/plume_regolith_solver_v1/src')

from scipy.integrate import quad
from scipy.optimize import fsolve

import var_timestepping as timestep
import var_range_of_interest as bounds
import var_soil as soil
import var_impinged_gas as imp
from fun_impinged_gas import *
from fun_soil import *
from fun_nozzle import *
from fun_integration import *


def compute_soil_density_depth(z):
    soil.rho_bulk = soil.rho_0 - (soil.rho_0 - soil.rho_inf)*(1-np.exp(-z/soil.h)) 
    return soil.rho_bulk

def integral_equation(h_new_exc, h_old_exc, mass_excavated):
    integral_equation = mass_excavated - quad(compute_soil_density_depth, h_old_exc, h_new_exc)[0]
    return integral_equation

h_old = 0
mass = 3279.544012872462

h_new = fsolve(integral_equation, h_old, args = (h_old, mass))
print(h_new)
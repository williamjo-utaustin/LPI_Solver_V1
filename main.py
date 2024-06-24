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
from fun_soil import *
from fun_nozzle import *
from fun_impinged_gas import *
from subroutines import *

# run main code
compute_nozzle_exhaust()
loop()

#plot_initial_disturbed_altitude()
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

print(bounds.r_centerlines_array_2)
plt.scatter(bounds.r_centerlines_array_2, np.ones(np.size(bounds.r_centerlines_array_2)))
plt.show()


#h_displaced_old[i] = h_displaced[i]
## assume it doubles and shit per timestep as a test 
#surf_den_excavated[i] = (m_excavated[i])/ring_area[i]
#
#def func(z):
#    y, err = quad(soil_density_depth, h_displaced[i], z)
#    return surf_den_excavated[i] - y
#
##print("Starting bound", h_displaced[i])
#sol = fsolve(func, 1E-8)
#h_displaced[i] = sol[0]
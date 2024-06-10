import numpy as np
import var_range_of_interest as bounds
from fun_soil import *

# universal variables
#rho_p = 1800
        
# Calculate the surface threshold parameters
E_th_0 = 0.123 #J/m^2/s
epsilon = 0.0029 # average value for soil efficiency
#alpha = 0.289 # J/m^3 (constant)
alpha_0 = 0.289
D_84 = 450E-6 # m, Check back with apollo 16 data for more accurate answer
D_bracket = 1.5 * D_84

rho_0 = 1100
rho_inf = 1800 #38% porosity deep within crater
h = 8.0/100 # From Metzger (2024)

k_cohesion = 50 #kg/m^3 (can change)

rho_bulk = None
E_downward = None

# compute the area of a differential ring
ring_area = np.zeros(bounds.n_points_centerline-1)
r_midpoint = np.zeros(bounds.n_points_centerline-1)

# compute the midpoint distances of each ring
m_dot_area_eroded = np.zeros(bounds.n_points_centerline)
m_dot_eroded_mid = np.zeros(bounds.n_points_centerline-1)
m_excavated_inst = np.zeros(bounds.n_points_centerline-1)
m_excavated_cumulative = np.zeros(bounds.n_points_centerline-1)

# initialize excavation heights
h_excavated_bounds = np.zeros(bounds.n_points_centerline)
h_excavated_mid = np.zeros(bounds.n_points_centerline-1)
h_excavated_mid_old = np.zeros(bounds.n_points_centerline-1)
surf_den_excavated = np.zeros(bounds.n_points_centerline-1)

# initialize densities
alpha = np.zeros(bounds.n_points_centerline)
E_th = np.zeros(bounds.n_points_centerline)

# initialize the midpoint value of the radial bins
for i in range(1,bounds.n_points_centerline):
    r_midpoint[i-1] = (bounds.r_centerlines_array[i] + bounds.r_centerlines_array[i-1])/2

# compute the ring areas
for i in range(1,bounds.n_points_centerline):
    ring_area[i-1] = np.pi * ((bounds.r_centerlines_array[i])**2 - (bounds.r_centerlines_array[i-1])**2) 
    #print(i-1, "to", i, bounds.r_centerlines_array[i-1], "to", bounds.r_centerlines_array[i], excavate.ring_area[i-1])

for i in range(0,bounds.n_points_centerline):
    alpha[i] = compute_cohesive_energy_density(0)
    E_th[i] = compute_threshold_energy(0)

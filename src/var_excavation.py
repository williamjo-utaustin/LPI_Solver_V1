#import numpy as np
#import var_range_of_interest as bounds
#from fun_soil import *
## compute the area of a differential ring
#ring_area = np.zeros(bounds.n_points_centerline-1)
#r_midpoint = np.zeros(bounds.n_points_centerline-1)
#
## compute the midpoint distances of each ring
#m_dot_area_eroded = np.zeros(bounds.n_points_centerline)
#m_dot_eroded_mid = np.zeros(bounds.n_points_centerline-1)
#m_excavated_inst = np.zeros(bounds.n_points_centerline-1)
#m_excavated_cumulative = np.zeros(bounds.n_points_centerline-1)
#
## initialize excavation heights
#h_excavated_bounds = np.zeros(bounds.n_points_centerline)
#h_excavated_mid = np.zeros(bounds.n_points_centerline-1)
#h_excavated_mid_old = np.zeros(bounds.n_points_centerline-1)
#surf_den_excavated = np.zeros(bounds.n_points_centerline-1)
#
## initialize densities
#alpha = np.zeros(bounds.n_points_centerline)
#E_th = np.zeros(bounds.n_points_centerline)
#
## initialize the midpoint value of the radial bins
#for i in range(1,bounds.n_points_centerline):
#    r_midpoint[i-1] = (bounds.r_centerlines_array[i] + bounds.r_centerlines_array[i-1])/2
#
## compute the ring areas
#for i in range(1,bounds.n_points_centerline):
#    ring_area[i-1] = np.pi * ((bounds.r_centerlines_array[i])**2 - (bounds.r_centerlines_array[i-1])**2) 
#    #print(i-1, "to", i, bounds.r_centerlines_array[i-1], "to", bounds.r_centerlines_array[i], excavate.ring_area[i-1])
#
#for i in range(0,bounds.n_points_centerline):
#    alpha[i] = compute_cohesive_energy_density(0)
#    E_th[i] = compute_threshold_energy(0)

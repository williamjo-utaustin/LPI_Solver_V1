import numpy as np
import sys

sys.path.append('/Users/williamjo/Documents/LPI/Codes/plume_regolith_solver_v1/src')

from fun_regolith_distributions import *
from fun_farfieldbounds import *
import var_soil as soil


plot_only = False


# determine the contribution of dust 
dust_wt_pct = compute_dust_mass_per_size_bin()

ejecta_ring_bounds =  np.genfromtxt("output/ejecta_ring_bounds.csv", delimiter = " ")

n_sim_particles_per_bin = 500
n_sorted_bins = 100
#range_of_analysis = 1000000
range_of_analysis = 1000

# index 1-9 [10], 1 - 9 microns
# index 10-19 [20], 10 - 90 microns
# index 19-29 [30], 100 - 900 microns
# index 30 - 36 [37], 1000 - 7000 microns

#index_grain_sizes = np.array([1,9])
#index_grain_sizes = np.array([10,18])
#index_grain_sizes = np.array([19,27])
#index_grain_sizes = np.array([28,36])

#index_grain_sizes = np.array([0,8])
index_grain_sizes = np.array([9,17])
#index_grain_sizes = np.array([18,26])
#index_grain_sizes = np.array([27,35])


area_sorted_bins, range_bounds, range_bounds_mid = create_far_field_bounds(n_sorted_bins, True)

#range_bounds = np.logspace(0, 7, n_sorted_bins+1)
#range_bounds = np.logspace(0, 4, n_sorted_bins+1)
#range_bounds = np.linspace(0, range_of_analysis, n_sorted_bins+1)
mass_upon_impact = np.zeros([np.size(soil.d_particle), n_sorted_bins])
mass_flux_upon_impact = np.zeros([np.size(soil.d_particle), n_sorted_bins])
momentum_flux_upon_impact = np.zeros([np.size(soil.d_particle), n_sorted_bins])
energy_flux_upon_impact = np.zeros([np.size(soil.d_particle), n_sorted_bins])
count_flux_upon_impact = np.zeros([np.size(soil.d_particle), n_sorted_bins])
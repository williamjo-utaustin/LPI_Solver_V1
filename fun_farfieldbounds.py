import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append('/Users/williamjo/Documents/LPI/Codes/plume_regolith_solver_v1/src')

import var_constants as cs

def magnitude(x):
    return np.sqrt(x[0]**2 + x[1]**2)

def create_far_field_bounds(n_segments, log_points):

    if(log_points):
        segment_bounds_y = np.logspace(-1, 6.54107976778, n_segments+1)
    else:
        segment_bounds_y = np.linspace(0, 2*cs.r_moon, n_segments+1)
    
    x_farfield_bounds = np.zeros_like(segment_bounds_y)
    y_farfield_bounds = np.zeros_like(segment_bounds_y)
    for i in range(0, np.size(segment_bounds_y)-1):

        y_farfield_bounds[i] = cs.r_moon - segment_bounds_y[i]
        x_farfield_bounds[i] = cs.r_moon * np.cos(np.arcsin(y_farfield_bounds[i]/cs.r_moon))

    x_farfield_bounds[np.size(segment_bounds_y)-1] = 0
    y_farfield_bounds[np.size(segment_bounds_y)-1] = -cs.r_moon


    ribbon_area = np.zeros(n_segments)
    for i in range(0, np.size(segment_bounds_y)-1):
        ribbon_height = np.abs(y_farfield_bounds[i+1] - y_farfield_bounds[i])
        ribbon_area[i] = 2*np.pi*cs.r_moon * ribbon_height

    particle_pos_init = np.array([x_farfield_bounds[0], y_farfield_bounds[0]])
    segment_bounds_arc = np.zeros_like(segment_bounds_y)
    for i in range(0, np.size(segment_bounds_y)-1):
        particle_pos = np.array([x_farfield_bounds[i], y_farfield_bounds[i]])
        alpha = np.arccos(np.dot(particle_pos_init, particle_pos)/(magnitude(particle_pos_init) * magnitude(particle_pos)))
        segment_bounds_arc[i] = alpha * cs.r_moon

    segment_bounds_arc[n_segments] = np.pi * cs.r_moon

    segment_midpoint_arc = np.zeros(n_segments)
    for i in range(0, n_segments):
        segment_midpoint_arc[i] = (segment_bounds_arc[i] + segment_bounds_arc[i+1])/2

    return ribbon_area, segment_bounds_arc, segment_midpoint_arc
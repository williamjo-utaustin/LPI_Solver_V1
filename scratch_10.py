import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append('/Users/williamjo/Documents/LPI/Codes/plume_regolith_solver_v1/src')

import var_constants as cs
#data_1 = np.genfromtxt("output/particle_range_0_1000.csv", delimiter= ",")
#data_2 = np.genfromtxt("output/particle_range_1010_2000.csv", delimiter= ",")
#data_3 = np.genfromtxt("output/particle_range_2050_2350.csv", delimiter= ",")
#data_4 = np.genfromtxt("output/particle_range_2360_2370.csv", delimiter= ",")
#def magnitude(x):
#    return np.sqrt(x[0]**2 + x[1]**2)
#
#
#def create_far_field_bounds(n_segments, log_points):
#
#    if(log_points):
#        segment_bounds_y = np.logspace(-1, 6.54107976778, n_segments+1)
#    else:
#        segment_bounds_y = np.linspace(0, 2*cs.r_moon, n_segments+1)
#    
#    x_farfield_bounds = np.zeros_like(segment_bounds_y)
#    y_farfield_bounds = np.zeros_like(segment_bounds_y)
#    for i in range(0, np.size(segment_bounds_y)-1):
#
#        y_farfield_bounds[i] = cs.r_moon - segment_bounds_y[i]
#        x_farfield_bounds[i] = cs.r_moon * np.cos(np.arcsin(y_farfield_bounds[i]/cs.r_moon))
#
#    x_farfield_bounds[np.size(segment_bounds_y)-1] = 0
#    y_farfield_bounds[np.size(segment_bounds_y)-1] = -cs.r_moon
#
#
#    ribbon_area = np.zeros(n_segments)
#    for i in range(0, np.size(segment_bounds_y)-1):
#        ribbon_height = np.abs(y_farfield_bounds[i+1] - y_farfield_bounds[i])
#        ribbon_area[i] = 2*np.pi*cs.r_moon * ribbon_height
#
#    particle_pos_init = np.array([x_farfield_bounds[0], y_farfield_bounds[0]])
#    segment_bounds_arc = np.zeros_like(segment_bounds_y)
#    for i in range(0, np.size(segment_bounds_y)-1):
#        particle_pos = np.array([x_farfield_bounds[i], y_farfield_bounds[i]])
#        alpha = np.arccos(np.dot(particle_pos_init, particle_pos)/(magnitude(particle_pos_init) * magnitude(particle_pos)))
#        segment_bounds_arc[i] = alpha * cs.r_moon
#
#    segment_bounds_arc[n_segments] = np.pi * cs.r_moon
#
#    segment_midpoint_arc = np.zeros(n_segments)
#    for i in range(0, n_segments):
#        segment_midpoint_arc[i] = (segment_bounds_arc[i] + segment_bounds_arc[i+1])/2
#
#    return ribbon_area, segment_bounds_arc, segment_midpoint_arc
#
#
#seg_area, seg_arc_bounds, seg_arc_mid = create_far_field_bounds(100, True)
#
#print(np.size(seg_area), np.size(seg_arc_bounds), np.size(seg_arc_mid))
#print(seg_area, seg_arc_bounds, seg_arc_mid)
#    
#
data_1 = np.genfromtxt('data/ejecta_ranges.csv', delimiter=',')




plt.semilogy(data_1[:,0], data_1[:,1])



#plt.semilogy(data_2[:,0], data_2[:,1])
#plt.semilogy(data_3[:,0], data_3[:,1])
#plt.semilogy(data_4[:,0], data_4[:,1])


sys.path.append('/Users/williamjo/Documents/LPI/Codes/plume_regolith_solver_v1/src')


import var_constants as cs

u_ej = np.linspace(0, 8000, 100)


v_ej_bar = (u_ej**2) / (cs.r_moon * cs.g)
r_max_1 = 2 * (cs.r_moon) * np.arctan((v_ej_bar**2 * np.sin(3*np.pi/180) * np.cos(3 * np.pi/180))/(1 - v_ej_bar**2 * (np.cos(3*np.pi/180))**2))


r_max_2 = (u_ej**2 * np.sin(2 * 3 * np.pi/180))/ (cs.g)

print(5442068.267824446*2)

for i in range(0,100):
    print(u_ej[i], r_max_1[i], r_max_2[i], 10920176)


vesc = np.sqrt(2 * cs.g * cs.r_moon)
print(vesc)
print(2 * np.pi * cs.r_moon)
plt.semilogy(u_ej, r_max_1)
plt.semilogy(u_ej, r_max_2)


print('{:.2e}'.format(1E-4 * (1/1800) * 1/((4/3)*np.pi * (0.5E-6)**3)))


plt.show()
import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg
import sys


sys.path.append('/Users/williamjo/Documents/LPI/Codes/plume_regolith_solver_v1/src')


from fun_regolith_distributions import *
import var_soil as soil
import var_solve_ejection as sol_ej
import var_constants as cs

# this will be the varying constant

def magnitude(x):
    return np.sqrt(x[0]**2 + x[1]**2)

def unit_norm(x):
    return x/np.sqrt(x[0]**2 + x[1]**2)

def g_vector(x):
    g_mag = cs.g * (cs.r_moon**2 / (x[0]**2 + x[1]**2))
    d_grav = -unit_norm(x)
    return d_grav * g_mag

# do not change these
delta_t = 1
launch_angle = 3*np.pi/180

min_vel = 1900
max_vel = 1901
v_init_list = np.linspace(min_vel,max_vel,1)
print(v_init_list)

save_values = np.zeros([1, 4])


for j in range(0,1):
    v_init = v_init_list[j]

    # set the initial values
    particle_pos = np.array([0, cs.r_moon])
    g_force = g_vector(particle_pos)
    particle_vel = np.array([v_init * np.cos(launch_angle), v_init * np.sin(launch_angle)])

    particle_pos_init = np.copy(particle_pos)

    # plot lunar circle
    x = np.linspace(-cs.r_moon,cs.r_moon, 1000)
    y = np.sqrt((cs.r_moon)**2 - x**2)

    plt.plot(x,y, color = 'black')
    plt.plot(x,-y, color = 'black')

    if(np.mod(j,1)==0):
        print(j, v_init)


    i = 0
    loop = True
    while loop:
        i = i + 1
        particle_vel = particle_vel + g_force * delta_t
        particle_pos = particle_pos + particle_vel * delta_t
        g_force = g_vector(particle_pos)


        #if(magnitude(particle_pos)-cs.r_moon > 8000):
        #    print(i * delta_t, particle_pos)
        #    break

        if(np.mod(i,10)==0):
            plt.scatter(particle_pos[0], particle_pos[1], s = 5)
            #print(particle_vel)
            print(i, magnitude(particle_pos) - cs.r_moon)
            #print(g_force)
        if(magnitude(particle_pos)<cs.r_moon or i > 100000000):
            #print("Iteration Total: ", i * delta_t)
            #print(particle_pos)
            #print(magnitude(particle_pos) - cs.r_moon)
            #print(particle_vel)

            alpha = np.arccos(np.dot(particle_pos_init, particle_pos)/(magnitude(particle_pos_init) * magnitude(particle_pos)))
            save_values[j, 0] = v_init
            save_values[j, 1] = alpha*cs.r_moon
            save_values[j, 2] = magnitude(particle_vel)
            save_values[j, 3] = i * delta_t
            
            #print(v_init, alpha * cs.r_moon, magnitude(particle_vel), i * delta_t)

            loop = False

#np.savetxt("output/particle_range_"+str(int(min_vel))+"_"+str(int(max_vel))+".csv", save_values, delimiter=",")

#plt.xlim(-2 * cs.r_moon, 2 * cs.r_moon)
#plt.ylim(-2 * cs.r_moon, 2 * cs.r_moon)
plt.gca().set_aspect('equal')
plt.show()

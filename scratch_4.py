import numpy as np
import matplotlib.pyplot as plt



def dudt(u_p, t, d_p, u_g, p_g, rho_g, T_g):
    
    A_const= 1.71575E-7
    beta = 0.78
    mfp = np.sqrt(np.pi/(2 * p_g * rho_g)) * A_const * T_g**beta
    
    Kn = mfp/d_p
    Re = (rho_g * np.abs(u_g - u_p) * d_p)/(A_const * T_g ** beta)
    
    if Re <= 0.0:
        CD = 0.0
    elif Re > 0.0 and Re <= 2:
        CD = (24.0/Re)
    elif Re >= 2 and Re < 500:
        CD = 18.5 * Re **-0.6
    else:
        CD = 0.44

    # add in correction factor to account for rarefied effects
    S_correction = 1 + Kn * (2.514 + 0.8 * np.exp(-0.55/Kn))
    CD = CD/S_correction

    C_constant = 3 * rho_g * CD / (4 * rho_p * d_p)
    #return C_constant * (u_g - u_p)**2
    return C_constant * (u_g - u_p)**2 + (1 - rho_g/rho_p) * 1.62

# adaptive timestepping
#r_d = 2400
#u_gas = 3665
#rho_p = 1100.35185
#rho_gas = 6.929E-11
#p_gas = 0.001376
#T_gas = 4762

#u_gas = 11.48
#rho_p = 1100.35185
#rho_gas = 1.49E-9
#p_gas = 0.0358
#T_gas = 4468

#r_d = 669
#u_gas = 4771.34
#rho_p = 1100.35185
#rho_gas = 1.13E-9
#p_gas = 464.389
#T_gas = 3706
#
#r_d = 1
#u_gas = 1607
#rho_p = 1100.35
#rho_gas = 3E-5
#p_gas = 716
#T_gas = 4423
#
#
#r_d = 10
#u_gas = 1591
#rho_p = 1100.35
#rho_gas = 1.3E-5
#p_gas = 100.38
#T_gas = 2586

# 1 m altitude test case
#r_d = 1
#u_gas = 3559
#rho_p = 1100.35
#rho_gas = 0.0008873
#p_gas = 11762
#T_gas = 3159.883
#data = np.genfromtxt('data/FontesEtAlDigitizedVelocity_40ton_1km.csv', delimiter= ',')


# 7 m altitude test case
r_d = 1
u_gas = 1300
rho_p = 1100.35
rho_gas = 3.01E-5
p_gas = 721
T_gas = 4439
data = np.genfromtxt('data/FontesEtAlDigitizedVelocity_40ton_7km.csv', delimiter= ',')

A_const= 1.71575E-7
beta = 0.78
#Re_g = rho_gas * u_gas * r_d /(A_const * T_gas**beta)
#print(Re_g)

#b_layer_thickness = 2.8 * r_d * (1/Re_g)**0.5
#print(b_layer_thickness)
d_particle_array = np.asarray([multiplier * magnitude for magnitude in [1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 1E-1] for multiplier in [1, 2, 3, 4, 5, 6, 7, 8, 9]])
#d_particle_array = np.asarray([multiplier * magnitude for magnitude in [1E-4, 1E-3, 1E-2, 1E-1] for multiplier in [1, 1.05, 1.10, 1.15, 1.20, 1.25, 1.5, 1.75, 1.9, 2, 5, 7, 9]])
u_ej = np.zeros_like(d_particle_array)
for i in range(0, np.size(d_particle_array)):
    x_p = 0
    y_p = 0
    
    u_p = 0
    t = 0

    # number of array points
    du = 1
    n_timesteps = (int(u_gas/du) - 1)

    t_array = np.zeros(n_timesteps)
    u_p_array = np.zeros(n_timesteps)
    
    x_p_array = np.zeros(n_timesteps)
    y_p_array = np.zeros(n_timesteps)
    
    d_p = d_particle_array[i]

    print(d_p) 
    for ts in range (0, n_timesteps):

        print(ts)
        u_p_array[ts] = u_p
        
        x_p_array[ts] = x_p
        y_p_array[ts] = y_p    
    
        t_array[ts] = t
        
        #if (y_p > 0.003):
        #    u_ej[i] = u_p
        #    break


        delta_t = du/(dudt(u_p, t, d_p, u_gas, p_gas, rho_gas, T_gas))

        u_p = u_p + du
        t = t + delta_t


        x_p = x_p + (u_p * delta_t) * np.cos(3 * np.pi/180)
        y_p = y_p + (u_p * delta_t) * np.sin(3 * np.pi/180)

        #print(t, x_p, y_p, u_p)
        plt.scatter(t, u_p, color = 'blue')

    plt.xlim(0,0.01)
    plt.ylim(0,1500)
    plt.xticks(fontsize = 18)
    plt.yticks(fontsize = 18)
    plt.xlabel("Time (s)", fontsize = 18)
    plt.ylabel("Velocity (m/s)", fontsize = 18)
    plt.show()
    


#plt.plot(d_p, x_p)
plt.loglog(data[:,0], data[:,1], linewidth = 3, label = 'Fontes et al. (2022)')
plt.loglog(d_particle_array, u_ej, linewidth = 3, label = 'ODE Computation')
plt.xlim(1E-6, 1E-3)
plt.xticks(fontsize = 18)
plt.yticks(fontsize = 18)
plt.xlabel("Particle Diameter (m)", fontsize = 18)
plt.ylabel("Ejecta Velocity (m/s)",fontsize = 18)
plt.legend(fontsize = 18)
#plt.ylim(1E0, 1E4)
plt.show()

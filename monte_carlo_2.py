import numpy as np
import matplotlib.pyplot as plt
import time
import sys

sys.path.append('/Users/williamjo/Documents/LPI/Codes/plume_regolith_solver_v1/src')



from fun_regolith_distributions import *
import var_soil as soil

plot_only = True



dust_wt_pct = compute_dust_mass_per_size_bin()

ejecta_ring_bounds =  np.genfromtxt("output/ejecta_ring_bounds.csv", delimiter = " ")

n_sim_particles_per_bin = 80
n_sorted_bins = 50
range_of_analysis = 200

# index 1-9 [10], 1 - 9 microns
# index 10-19 [20], 10 - 90 microns
# index 19-29 [30], 100 - 900 microns
# index 30 - 36 [37], 1000 - 7000 microns

#index_grain_sizes = np.array([0,9])
#index_grain_sizes = np.array([9,18])
#index_grain_sizes = np.array([18,27])
index_grain_sizes = np.array([27,36])

range_bounds = np.logspace(0, 2, n_sorted_bins+1)
#range_bounds = np.linspace(0, range_of_analysis, n_sorted_bins+1)
mass_upon_impact = np.zeros([np.size(soil.d_particle), n_sorted_bins])
mass_flux_upon_impact = np.zeros([np.size(soil.d_particle), n_sorted_bins])
momentum_flux_upon_impact = np.zeros([np.size(soil.d_particle), n_sorted_bins])
energy_flux_upon_impact = np.zeros([np.size(soil.d_particle), n_sorted_bins])
area_sorted_bins = np.zeros(n_sorted_bins)
range_bounds_mid = np.zeros(n_sorted_bins)

for i in range(0, n_sorted_bins):
    area_sorted_bins[i] = np.pi * (range_bounds[i+1]**2 - range_bounds[i]**2)
    range_bounds_mid[i] = (range_bounds[i+1] + range_bounds[i])/2



if(plot_only is False):

    np.savetxt("output/range_bounds_"+str(index_grain_sizes[0])+"_"+str(index_grain_sizes[1])+".csv", range_bounds_mid, delimiter=",")

    # -----------------------------------------------            
    # loop through particle launch time
    # -----------------------------------------------            
    for t in range(0,301):

        if(np.mod(t,100)==0): 
            print("Time (s) = ", t * 0.01)

        ejecta_properties = np.genfromtxt("output/ejecta_properties_"+str(t)+".csv", delimiter = " ")
        ejecta_velocities = np.genfromtxt("output/ejecta_velocities_"+str(t)+".csv", delimiter = " ")


        # extract the instantaneous mass excavated data for each radial ring at time i * dt
        ejecta_indices = ejecta_properties[:,0]
        ejecta_mass_excavated = ejecta_properties[:,3]/1# Set to 0.25 degree ring
        ejecta_rho_excavated = ejecta_properties[:,2] # Set to 0.25 degree ring

        # set 10 cells per size bin spaced evenly apart
        u_ej_sim_particle_mid = np.zeros(np.size(ejecta_velocities[:,1])+2)
        r_ej_sim_particle_mid = np.zeros(np.size(ejecta_velocities[:,1])+2)
        r_ej_sim_particle_bounds_lower = np.array(ejecta_ring_bounds[0:np.size(ejecta_velocities[:,1]), 1])
        r_ej_sim_particle_bounds_upper = np.array(ejecta_ring_bounds[0:np.size(ejecta_velocities[:,1]), 3])

        #for i in range(0, np.size(ejecta_velocities[:,1])):
        #    print(i, r_ej_sim_particle_bounds_lower[i], r_ej_sim_particle_bounds_upper[i])



        # ---------------------------------------------------------- 
        # loop through all particle sizes
        # ---------------------------------------------------------- 
        for i in range(index_grain_sizes[0],index_grain_sizes[1]):   

            proportion = dust_wt_pct[i-1]

            #print("Particle Size", soil.d_particle[i-1])
            #print("Proportion", proportion)

            # --------------------------------------------------------------- 
            # extract the ring bounds and the ejecta velocities from datasets
            # --------------------------------------------------------------- 
            for j in range(1,np.size(u_ej_sim_particle_mid)-1):
                r_ej_sim_particle_mid[j] = ejecta_ring_bounds[j-1,2]
                u_ej_sim_particle_mid[j] = ejecta_velocities[j-1,i]
                #print(r_ej_sim_particle_mid[j])
                #print(u_ej_sim_particle_mid[j])


            r_ej_sim_particle_mid[np.size(u_ej_sim_particle_mid)-1] = ejecta_ring_bounds[np.size(u_ej_sim_particle_mid)-3,3]
            #print(r_ej_sim_particle_mid, np.size(r_ej_sim_particle_mid))
            #print(u_ej_sim_particle_mid, np.size(u_ej_sim_particle_mid))


            u_ej_sim_particle_bounds_lower = np.interp(r_ej_sim_particle_bounds_lower, r_ej_sim_particle_mid, u_ej_sim_particle_mid)
            u_ej_sim_particle_bounds_upper = np.interp(r_ej_sim_particle_bounds_upper, r_ej_sim_particle_mid, u_ej_sim_particle_mid)


            # ---------------------------------------------------------- 
            # loop through each launch bin range
            # ---------------------------------------------------------- 
            for j in range(0, np.size(u_ej_sim_particle_bounds_lower)):
            
                #print(j, r_ej_sim_particle_bounds_lower[j], r_ej_sim_particle_bounds_upper[j], u_ej_sim_particle_bounds_lower[j], u_ej_sim_particle_bounds_upper[j], ejecta_mass_excavated[j])

                r_ej_sim_particle = np.linspace(r_ej_sim_particle_bounds_lower[j], r_ej_sim_particle_bounds_upper[j], n_sim_particles_per_bin)
                u_ej_sim_particle = np.linspace(u_ej_sim_particle_bounds_lower[j], u_ej_sim_particle_bounds_upper[j], n_sim_particles_per_bin)


                total_mass_excavated_particle_size = ejecta_mass_excavated[j] * dust_wt_pct[i-1]

                #print(" ") 
                #print("Launch Bin #", j)
                #print("Total Mass Excavated Within this Range (kg)", ejecta_mass_excavated[j])
                #print("Mass Excavated for Particle Size ", soil.d_particle[i-1], "m (kg): ", total_mass_excavated_particle_size)
                #print(" ")
                #print("Distance from Centerline (m), Speed of Launch (m/s), Landing Location (m), sim particle mass (kg)")

                if(j == 0):
                    lower_bound = 1
                    upper_bound = np.size(r_ej_sim_particle)
                    sim_particles = n_sim_particles_per_bin - 1  

                elif (j == np.size(u_ej_sim_particle_bounds_lower)-1):
                    lower_bound = 0
                    upper_bound = np.size(r_ej_sim_particle) - 1
                    sim_particles = n_sim_particles_per_bin - 1  

                else:
                    lower_bound = 0
                    upper_bound = np.size(r_ej_sim_particle)
                    sim_particles = n_sim_particles_per_bin

            # ---------------------------------------------------------- 
            # loop through each sim particle in each bin
            # ---------------------------------------------------------- 
                for k in range(lower_bound, upper_bound):


                    #print(u_ej_sim_particle[k])
                    #v_ej_bar = u_ej_sim_particle[k]**2 / (1738E3 * 1.62)
                    #print("v_ej", v_ej_bar)
                    #r_max_2 = 2 * (1738E3) * np.arctan((v_ej_bar**2 * np.sin(3*np.pi/180) * np.cos(3 * np.pi/180))/(1-v_ej_bar**2 * (np.cos(3*np.pi/180))**2))


                    r_max = (u_ej_sim_particle[k]**2 * np.sin(2 * 3 * np.pi/180))/ (1.62)
                    #print("Rmax", r_max, r_max_2)

                    mass_sim_particle = total_mass_excavated_particle_size/sim_particles

                    # sort the sim particle masses into bins based on landing location
                    for l in range(1, n_sorted_bins):

                        if(range_bounds[l] > r_max + r_ej_sim_particle[k]):


                            #if (l > 11): 
                            #    print(k, r_ej_sim_particle[k],  r_max + r_ej_sim_particle[k], l-1, u_ej_sim_particle[k], mass_sim_particle)


                            mass_upon_impact[i-1,l-1] = mass_upon_impact[i-1,l-1] + mass_sim_particle
                            momentum_flux_upon_impact[i-1,l-1] = momentum_flux_upon_impact[i-1,l-1] + (mass_sim_particle * u_ej_sim_particle[k])/area_sorted_bins[l-1]
                            mass_flux_upon_impact[i-1,l-1] = mass_flux_upon_impact[i-1,l-1] + (mass_sim_particle/area_sorted_bins[l-1])
                            energy_flux_upon_impact[i-1,l-1] = energy_flux_upon_impact[i-1,l-1] + (0.5 * mass_sim_particle * u_ej_sim_particle[k]**2)/area_sorted_bins[l-1]
                            break
                        
                #print(mass_flux_upon_impact[i-1,:])
                #print(energy_flux_upon_impact[i-1,:])
                #print(energy_flux_upon_impact[i-1,:])


    np.savetxt("output/mass_flux_upon_impact_"+str(index_grain_sizes[0])+"_"+str(index_grain_sizes[1])+".csv", mass_flux_upon_impact, delimiter=",")
    np.savetxt("output/momentum_flux_upon_impact_"+str(index_grain_sizes[0])+"_"+str(index_grain_sizes[1])+".csv", momentum_flux_upon_impact, delimiter=",")
    np.savetxt("output/energy_flux_upon_impact_"+str(index_grain_sizes[0])+"_"+str(index_grain_sizes[1])+".csv", energy_flux_upon_impact, delimiter=",")


if(plot_only):
    range_bounds_mid = np.genfromtxt("output/range_bounds_"+str(index_grain_sizes[0])+"_"+str(index_grain_sizes[1])+".csv", delimiter=",")
    mass_flux_upon_impact = np.genfromtxt("output/mass_flux_upon_impact_"+str(index_grain_sizes[0])+"_"+str(index_grain_sizes[1])+".csv", delimiter=',')
    momentum_flux_upon_impact = np.genfromtxt("output/momentum_flux_upon_impact_"+str(index_grain_sizes[0])+"_"+str(index_grain_sizes[1])+".csv", delimiter=',')
    energy_flux_upon_impact = np.genfromtxt("output/energy_flux_upon_impact_"+str(index_grain_sizes[0])+"_"+str(index_grain_sizes[1])+".csv", delimiter=',')

total_mass_flux_upon_impact = np.zeros(n_sorted_bins)
total_momentum_flux_upon_impact = np.zeros(n_sorted_bins)
total_energy_flux_upon_impact = np.zeros(n_sorted_bins)

for i in range(0,n_sorted_bins):
    total_mass_flux_upon_impact[i] = np.sum(mass_flux_upon_impact[:,i])
    total_momentum_flux_upon_impact[i] = np.sum(momentum_flux_upon_impact[:,i])
    total_energy_flux_upon_impact[i] = np.sum(energy_flux_upon_impact[:,i])

plt.figure(figsize=(5,5))
plt.plot(range_bounds_mid, total_energy_flux_upon_impact, color = 'black', linewidth = 2, label = 'Total')

for i in range(index_grain_sizes[0],index_grain_sizes[1]):
    #plt.ylim(0,10)
    plt.plot(range_bounds_mid, energy_flux_upon_impact[i,:], linewidth = 2, linestyle = 'dotted', label = "Particle Size (m): " + str('{:.1E}'.format(soil.d_particle[i])))
    #plt.xscale('log')
    plt.yscale('log')

plt.xlabel("Distance from Centerline (m)", fontsize = 16)
plt.ylabel("Impact Energy per Unit Area $(J/m^2)$", fontsize = 16)
plt.xticks(fontsize = 18)
plt.yticks(fontsize = 18)
plt.locator_params(axis='x', nbins=5)
plt.xlim(0,range_of_analysis)
#plt.ylim(1E-1,1E3)
plt.grid()
plt.legend(fontsize = 8)
plt.savefig("output/impact_energy_profile_"+str(index_grain_sizes[0])+"_"+str(index_grain_sizes[1])+".png", bbox_inches='tight', dpi=100)
plt.close()


plt.figure(figsize=(5,5))
plt.plot(range_bounds_mid, total_momentum_flux_upon_impact, color = 'black', linewidth = 2, label = 'Total')

for i in range(index_grain_sizes[0], index_grain_sizes[1]):
    #plt.ylim(0,10)
    plt.plot(range_bounds_mid, momentum_flux_upon_impact[i,:], linewidth = 2, linestyle = 'dotted', label = "Particle Size (m): " + str('{:.1E}'.format(soil.d_particle[i])))
    #plt.xscale('log')
    plt.yscale('log')

plt.xlabel("Distance from Centerline (m)", fontsize = 16)
plt.ylabel("Impact Momentum per Unit Area $(N \cdot s/m^2)$", fontsize = 16)
plt.xticks(fontsize = 18)
plt.yticks(fontsize = 18)
plt.locator_params(axis='x', nbins=5)
plt.xlim(0,range_of_analysis)
plt.grid()
plt.legend(fontsize = 8)
plt.savefig("output/impact_momentum_profile_"+str(index_grain_sizes[0])+"_"+str(index_grain_sizes[1])+".png", bbox_inches='tight', dpi=100)
plt.close()

plt.figure(figsize=(5,5))
plt.plot(range_bounds_mid, total_mass_flux_upon_impact, color = 'black', linewidth = 2, label = 'Total')

for i in range(index_grain_sizes[0],index_grain_sizes[1]):
    #plt.ylim(0,10)
    plt.plot(range_bounds_mid, mass_flux_upon_impact[i,:], linewidth = 2, linestyle = 'dotted', label = "Particle Size (m): " + str('{:.1E}'.format(soil.d_particle[i])))
    #plt.xscale('log')
    plt.yscale('log')

plt.xlabel("Distance from Centerline (m)", fontsize = 16)
plt.ylabel("Impact Mass per Unit Area $(kg/m^2)$", fontsize = 16)
plt.xticks(fontsize = 18)
plt.yticks(fontsize = 18)
plt.locator_params(axis='x', nbins=5)
plt.xlim(0,range_of_analysis)
plt.grid()
plt.legend(fontsize = 8)
plt.savefig("output/impact_mass_profile_"+str(index_grain_sizes[0])+"_"+str(index_grain_sizes[1])+".png", bbox_inches='tight',dpi=100)
plt.close()







#print(u_ej_sim_particle_bounds)

        #for j in range(1, np.size(u_ej_sim_particle_bounds)):
        #    #print(j-1, u_ej_sim_particle_bounds[j-1], u_ej_sim_particle_bounds[j])
        #    r_ej_sim_particle = np.linspace(r_ej_sim_particle_bounds[j-1], r_ej_sim_particle_bounds[j], int(n_sim_particles_per_bin))
        #    u_ej_sim_particle = np.linspace(u_ej_sim_particle_bounds[j-1], u_ej_sim_particle_bounds[j], int(n_sim_particles_per_bin))

        #    print("Set ", j) 
        #    for k in range(0, np.size(r_ej_sim_particle)):
        #        r_max = (u_ej_sim_particle[k]**2 * np.sin(2 * 3 * np.pi/180))/ (1.62)
        #        print(r_ej_sim_particle[k], u_ej_sim_particle[k], r_max + r_ej_sim_particle[k])
        #    print(" ")
    
            #dust_mass, n_dust_particles_total = compute_dust_mass_per_size_bin(ejecta_mass_excavated[j-1],ejecta_rho_excavated[j-1])







            #print(r_ej_sim_particle)



        #print((r_ej_sim_particle_bounds)) 
        #print((u_ej_sim_particle_bounds))

        #plt.scatter(r_ej_sim_particle_bounds, u_ej_sim_particle_bounds)
        #plt.show()











    # maybe we set 10 simulation particles per ring
    #for j in range(0, np.size(ejecta_mass_excavated)+1):

    #    print(j)

    #    if(j == 0): 
    #        r_sim_particle = np.linspace(ejecta_ring_bounds[j,1], ejecta_ring_bounds[j,2], int(n_sim_particles_per_bin/2))
    #        u_ej_sim_particle = np.interp(r_sim_particle, np.array([ejecta_ring_bounds[j,1], ejecta_ring_bounds[j,2]]), np.array([0, ejecta_velocities[j,1]]))
    #        for k in range(1, np.size(u_ej_sim_particle)): 
    #            r_max = (u_ej_sim_particle[k]**2 * np.sin(3 * np.pi/180))/ (1.62)
    #            print(r_sim_particle[k], u_ej_sim_particle[k], r_max + r_sim_particle[k])
    #        
    #    
    #    elif(j == np.size(ejecta_mass_excavated)):
    #        r_sim_particle = np.linspace(ejecta_ring_bounds[j-1,2], ejecta_ring_bounds[j-1,3], int(n_sim_particles_per_bin/2))
    #        u_ej_sim_particle = np.interp(r_sim_particle, np.array([ejecta_ring_bounds[j-1,2], ejecta_ring_bounds[j-1,3]]), np.array([ejecta_velocities[np.size(ejecta_mass_excavated)-1,1], 0]))
    #        for k in range(0, np.size(u_ej_sim_particle)-1): 
    #            r_max = (u_ej_sim_particle[k]**2 * np.sin(3 * np.pi/180))/ (1.62)
    #            print(r_sim_particle[k], u_ej_sim_particle[k], r_max + r_sim_particle[k])

    #    else:
    #        r_sim_particle = np.linspace(ejecta_ring_bounds[j-1,2], ejecta_ring_bounds[j,2], int(n_sim_particles_per_bin))
    #        u_ej_sim_particle = np.interp(r_sim_particle, np.array([ejecta_ring_bounds[j-1,2], ejecta_ring_bounds[j,2]]), np.array([ejecta_velocities[j-1,1], ejecta_velocities[j,1]]))

    #        for k in range(0, np.size(u_ej_sim_particle)): 
    #            r_max = (u_ej_sim_particle[k]**2 * np.sin(3 * np.pi/180))/ (1.62)
    #            print(r_sim_particle[k], u_ej_sim_particle[k], r_max + r_sim_particle[k])
    #        
    #        
    #        
    #        
    #        #print(r_sim_particle)
    #    print("   ")
    #    #for k in range(0, np.size(r_sim_particle)):




    #    #print(ejecta_ring_bounds[j,1], ejecta_ring_bounds[j,2], ejecta_ring_bounds[j,3], ejecta_velocities[j,1], r_max + ejecta_ring_bounds[j,2])












    



    ## loop through each ring

    #for j in range(0, np.size(ejecta_mass_excavated)):
    ##for j in range(0, 4):
    #    
    #    
    #    dust_mass, n_dust_particles_total = compute_dust_mass_per_size_bin(ejecta_mass_excavated[j],ejecta_rho_excavated[j])
    #    
    #    sim_dust_particles = np.ones(np.size(dust_mass)) # we set this to 100 sim particles per size bin but we can change this later!
    #    for k in range(0, np.size(dust_mass)):
    #        sim_dust_particles[k] = sim_dust_particles[k] * int(250 *(1/(k+1)))

    #    f_num_grains = n_dust_particles_total/sim_dust_particles
    #    
    #    
    #    print("                ") 
    #    print("ring count", j)
    #    print("Mass Excavated in ring j", ejecta_mass_excavated[j])
    #    print("simulation particles per grain size", sim_dust_particles)
    #    print("total dust particles", np.sum(sim_dust_particles)) 
    #    print("f_num", f_num_grains)



    #    # begin launching particles for each ring
    #    
    #    # set launch positions
    #    r_p = np.zeros([int(sum(sim_dust_particles)), 3])

    #    # set launch velocities
    #    v_p = np.zeros([int(sum(sim_dust_particles)), 3]) 

    #    index_count = 0
    #    size_bin = 0

import numpy as np
import matplotlib.pyplot as plt
import time
import sys

sys.path.append('/Users/williamjo/Documents/LPI/Codes/plume_regolith_solver_v1/src')

import var_soil as soil


np.random.seed(67)

def compute_dust_mass_per_size_bin():


    data = np.genfromtxt('/Users/williamjo/Documents/LPI/Raw Data/a16_dataset.csv', delimiter=',') * 1E-6
    #print(data)

    #plt.semilogx(data[:,0],data[:,1])
    #plt.semilogx(data[:,2],data[:,3])

    d_set_1_y = np.empty(np.size(soil.d_particle))
    d_set_2_y = np.empty(np.size(soil.d_particle))

    d_set_1_y = np.interp(soil.d_particle, data[:,0], data[:,1])
    d_set_2_y = np.interp(soil.d_particle, data[:,2], data[:,3])
    d_set_3_y = (d_set_1_y + d_set_2_y)/2

    d_set_3_y = d_set_3_y*100/np.max(d_set_3_y)

    # create the weight percent for each 
    wt_pct = np.empty(np.size(soil.d_particle))
    count = 0

    # print out the index, radius range (0 to 1 microns) for a bin, and then the wt% of this range
    wt_pct_cumulative =  d_set_3_y[0]
    
    #print(count, 0,"to", d_particle[0], d_set_3_y[0])
    sum_wt = d_set_3_y[0]
    wt_pct[0] = d_set_3_y[0]

    # print out the index, radius range (1 to [900,950, or < 1000] microns) for a bin, and then the wt% of this range
    for i in range(1,np.size(soil.d_particle)):
        count = count + 1
        sum_wt = sum_wt + d_set_3_y[i] - d_set_3_y[i-1]
        wt_pct[i] = d_set_3_y[i] - d_set_3_y[i-1]
        wt_pct_cumulative =  wt_pct_cumulative + wt_pct[i]
        #print(count, d_particle[i-1],"to",d_particle[i],d_set_3_y[i] - d_set_3_y[i-1], wt_pct_cumulative)

    # convert the weight percent to decimals (sum should equal 1)
    wt_pct = wt_pct/100
    #print(wt_pct)
    #print(wt_pct.size)

    # given the particle mass, use the distributed wt% of the ranges to obtain the total mass of the grains in each size bin
    # note that dust_mass will change further down the code (in this case, it represents the original ice )
    #dust_mass = m_parcel * wt_pct

    # compute the mass of each individual particle sizes
    #mass_per_particle_set_1_x = avg_den_excavated * (4/3) * np.pi * (d_particle/2)**3
    #total_dust_particles_per_bin = dust_mass/mass_per_particle_set_1_x


    #return dust_mass, total_dust_particles_per_bin

    return wt_pct


dust_wt_pct = compute_dust_mass_per_size_bin()

print(dust_wt_pct, np.size(dust_wt_pct))


ejecta_ring_bounds =  np.genfromtxt("output/ejecta_ring_bounds.csv", delimiter = " ")

n_sim_particles_per_bin = 100


n_sorted_bins = 200
#range_bounds = np.logspace(0, 2, n_sorted_bins+1)
range_bounds = np.linspace(0, 10000, n_sorted_bins+1)
mass_upon_impact = np.zeros([np.size(soil.d_particle), n_sorted_bins])
mass_flux_upon_impact = np.zeros([np.size(soil.d_particle), n_sorted_bins])
energy_flux_upon_impact = np.zeros([np.size(soil.d_particle), n_sorted_bins])
area_sorted_bins = np.zeros(n_sorted_bins)
range_bounds_mid = np.zeros(n_sorted_bins)

for i in range(0, n_sorted_bins):
    area_sorted_bins[i] = np.pi * (range_bounds[i+1]**2 - range_bounds[i]**2)
    range_bounds_mid[i] = (range_bounds[i+1] + range_bounds[i])/2

# -----------------------------------------------            
# loop through particle launch time
# -----------------------------------------------            
for t in range(0,301):
    
    print("Timestep", t)
    
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
    for i in range(1,2):   
        
        proportion = dust_wt_pct[i-1]

        print("Particle Size", soil.d_particle[i-1])
        print("Proportion", proportion)

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

            print(" ") 
            print("Launch Bin #", j)
            print("Total Mass Excavated Within this Range (kg)", ejecta_mass_excavated[j])
            print("Mass Excavated for Particle Size ", soil.d_particle[i-1], "m (kg): ", total_mass_excavated_particle_size)
            print(" ")
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
                v_ej_bar = u_ej_sim_particle[k]**2 / (1738E3 * 1.62)
                #print("v_ej", v_ej_bar)
                #r_max = 2 * (1738E3) * np.arctan((v_ej_bar**2 * np.sin(3*np.pi/180) * np.cos(3 * np.pi/180))/(1-v_ej_bar**2 * (np.cos(3*np.pi/180))**2))
                r_max = (u_ej_sim_particle[k]**2 * np.sin(2 * 3 * np.pi/180))/ (1.62)
                print("Rmax", r_max)
                mass_sim_particle = total_mass_excavated_particle_size/sim_particles
            
                # sort the sim particle masses into bins based on landing location
                for l in range(1, n_sorted_bins):
                    
                    if(range_bounds[l] > r_max + r_ej_sim_particle[k]):


                        if (l > 11): 
                            print(k, r_ej_sim_particle[k],  r_max + r_ej_sim_particle[k], l-1, u_ej_sim_particle[k], mass_sim_particle)
                        
                        
                        mass_upon_impact[i-1,l-1] = mass_upon_impact[i-1,l-1] + mass_sim_particle
                        mass_flux_upon_impact[i-1,l-1] = mass_flux_upon_impact[i-1,l-1] + (mass_sim_particle/area_sorted_bins[l-1])
                        energy_flux_upon_impact[i-1,l-1] = energy_flux_upon_impact[i-1,l-1] + (0.5 * mass_sim_particle * u_ej_sim_particle[k]**2)/area_sorted_bins[l-1]
                        break
            
            #print(mass_flux_upon_impact[i-1,:])
            print(energy_flux_upon_impact[i-1,:])
            #print(energy_flux_upon_impact[i-1,:])
    
    if(np.mod(t,100)==0):
        plt.title("Time = "+str(t))
        #plt.ylim(0,10)
        plt.scatter(range_bounds_mid, energy_flux_upon_impact[i-1,:])
        #plt.xscale('log')
        #plt.yscale('log')
        plt.show()
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

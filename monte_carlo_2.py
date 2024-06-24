import numpy as np
import matplotlib.pyplot as plt
import time 

np.random.seed(67)

def compute_dust_mass_per_size_bin(m_parcel, avg_den_excavated):

    # particles/bin
    comp_particles_per_size_bin = 100

    data = np.genfromtxt('/Users/williamjo/Documents/LPI/Raw Data/a16_dataset.csv', delimiter=',') * 1E-6
    #print(data)

    #plt.semilogx(data[:,0],data[:,1])
    #plt.semilogx(data[:,2],data[:,3])

    d_particle = np.asarray([multiplier * magnitude for magnitude in [1E-6, 1E-5, 1E-4, 1E-3] for multiplier in [1, 2, 3, 4, 5, 6, 7, 8, 9]])
    d_set_1_y = np.empty(d_particle.size)
    d_set_2_y = np.empty(d_particle.size)

    d_set_1_y = np.interp(d_particle, data[:,0], data[:,1])
    d_set_2_y = np.interp(d_particle, data[:,2], data[:,3])
    d_set_3_y = (d_set_1_y + d_set_2_y)/2

    d_set_3_y = d_set_3_y*100/np.max(d_set_3_y)

    # create the weight percent for each 
    wt_pct = np.empty(d_particle.size)
    count = 0

    # print out the index, radius range (0 to 1 microns) for a bin, and then the wt% of this range
    wt_pct_cumulative =  d_set_3_y[0]
    
    #print(count, 0,"to", d_particle[0], d_set_3_y[0])
    sum_wt = d_set_3_y[0]
    wt_pct[0] = d_set_3_y[0]

    # print out the index, radius range (1 to [900,950, or < 1000] microns) for a bin, and then the wt% of this range
    for i in range(1, d_particle.size):
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
    dust_mass = m_parcel * wt_pct

    # compute the mass of each individual particle sizes
    mass_per_particle_set_1_x = avg_den_excavated * (4/3) * np.pi * (d_particle/2)**3
    total_dust_particles_per_bin = dust_mass/mass_per_particle_set_1_x


    return dust_mass, total_dust_particles_per_bin

ejecta_ring_bounds =  np.genfromtxt("output/ejecta_ring_bounds.csv", delimiter = " ")

n_sim_particles_per_bin = 10
# loop through particle launch time
for t in range(0,1):
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
    
    for i in range(0, np.size(ejecta_velocities[:,1])):
        print(i, r_ej_sim_particle_bounds_lower[i], r_ej_sim_particle_bounds_upper[i])
    

    for i in range(1,2):   
    
    #for i in range(1,37):


        for j in range(1,np.size(u_ej_sim_particle_mid)-1):
            r_ej_sim_particle_mid[j] = ejecta_ring_bounds[j-1,2]
            u_ej_sim_particle_mid[j] = ejecta_velocities[j-1,i]
        

        r_ej_sim_particle_mid[np.size(u_ej_sim_particle_mid)-1] = ejecta_ring_bounds[np.size(u_ej_sim_particle_mid)-3,3]
        #print(r_ej_sim_particle_mid, np.size(r_ej_sim_particle_mid))
        #print(u_ej_sim_particle_mid, np.size(u_ej_sim_particle_mid))


        u_ej_sim_particle_bounds_lower = np.interp(r_ej_sim_particle_bounds_lower, r_ej_sim_particle_mid, u_ej_sim_particle_mid)
        u_ej_sim_particle_bounds_upper = np.interp(r_ej_sim_particle_bounds_upper, r_ej_sim_particle_mid, u_ej_sim_particle_mid)
        
        for j in range(0, np.size(u_ej_sim_particle_bounds_lower)):
        
            #print(j, r_ej_sim_particle_bounds_lower[j], r_ej_sim_particle_bounds_upper[j], u_ej_sim_particle_bounds_lower[j], u_ej_sim_particle_bounds_upper[j], ejecta_mass_excavated[j])

            r_ej_sim_particle = np.linspace(r_ej_sim_particle_bounds_lower[j], r_ej_sim_particle_bounds_upper[j], int(n_sim_particles_per_bin))
            u_ej_sim_particle = np.linspace(u_ej_sim_particle_bounds_lower[j], u_ej_sim_particle_bounds_upper[j], int(n_sim_particles_per_bin))


            print("Launch Bin #", j)
            print("Total Mass Excavated Within this Range (kg)", ejecta_mass_excavated[j])
            print("Mass Excavated for Particle Size: ", i)
            print("Distance from Centerline (m), Speed of Launch (m/s), Landing Location (m)")



            for k in range(0, np.size(r_ej_sim_particle)):
                r_max = (u_ej_sim_particle[k]**2 * np.sin(2 * 3 * np.pi/180))/ (1.62)
                print(r_ej_sim_particle[k], u_ej_sim_particle[k], r_max + r_ej_sim_particle[k])
            print(" ")




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







    time.sleep(100)





    



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

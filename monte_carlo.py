import numpy as np
import matplotlib.pyplot as plt

np.random.seed(67)


n_particles = 1000


r_inner = 100
r_outer = 105

r_particles = np.random.uniform(r_inner, r_outer, n_particles)
phi_launch = np.random.uniform(0, np.pi/4, n_particles)

#print(r_particles, phi_launch)

x = r_particles * np.cos(phi_launch)
y = r_particles * np.sin(phi_launch)


#plt.scatter(x,y)
#plt.xlim(-120,120)
#plt.ylim(-120,120)
#plt.show()

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



# simulate the ballistics in 0.25 degree slice

ejecta_ring_bounds =  np.genfromtxt("output/ejecta_ring_bounds.csv", delimiter = " ")


# loop through particle launch time
for i in range(0,1):
    ejecta_properties = np.genfromtxt("output/ejecta_properties_"+str(i)+".csv", delimiter = " ")
    ejecta_velocities = np.genfromtxt("output/ejecta_velocities_"+str(i)+".csv", delimiter = " ")

    # extract the instantaneous mass excavated data for each radial ring at time i * dt
    ejecta_indices = ejecta_properties[:,0]
    ejecta_mass_excavated = ejecta_properties[:,3]/1# Set to 0.25 degree ring
    ejecta_rho_excavated = ejecta_properties[:,2] # Set to 0.25 degree ring

    # loop through each ring

    #for j in range(0, np.size(ejecta_mass_excavated)-1):
    for j in range(0, 4):
        
        
        dust_mass, n_dust_particles_total = compute_dust_mass_per_size_bin(ejecta_mass_excavated[j],ejecta_rho_excavated[j])
        
        sim_dust_particles = np.ones(np.size(dust_mass)) # we set this to 100 sim particles per size bin but we can change this later!
        for k in range(0, np.size(dust_mass)):
            sim_dust_particles[k] = sim_dust_particles[k] * int(250 *(1/(k+1)))

        f_num_grains = n_dust_particles_total/sim_dust_particles
        
        
        #print("                ") 
        #print("ring count", j)
        #print("Mass Excavated in ring j", ejecta_mass_excavated[j])
        #print("simulation particles per grain size", sim_dust_particles)
        #print("total dust particles", np.sum(sim_dust_particles)) 
        #print("f_num", f_num_grains)



        # begin launching particles for each ring
        
        # set launch positions
        r_p = np.zeros([int(sum(sim_dust_particles)), 3])

        # set launch velocities
        v_p = np.zeros([int(sum(sim_dust_particles)), 3]) 

        index_count = 0
        size_bin = 0



        for k in range(0,int(sum(sim_dust_particles))):

            # set launch positions
            r_particles = np.random.uniform(ejecta_ring_bounds[int(ejecta_indices[j]),1], ejecta_ring_bounds[int(ejecta_indices[j]),3]) 
            phi_launch = np.random.uniform(0, 2*np.pi/1440)

            # launch position xp
            r_p[k,0] = r_particles * np.cos(phi_launch)
            r_p[k,1] = r_particles * np.sin(phi_launch)
            r_p[k,2] = 0.3

            # launch velocity vp
            v_p[k,0] = ejecta_velocities[j,size_bin+1] * np.cos(3 * np.pi/180) * np.cos(phi_launch)
            v_p[k,1] = ejecta_velocities[j,size_bin+1] * np.cos(3 * np.pi/180) * np.sin(phi_launch)
            v_p[k,2] = ejecta_velocities[j,size_bin+1] * np.sin(3 * np.pi/180)

            print(k, ejecta_velocities[j,size_bin], v_p[k,0], v_p[k,1], v_p[k,2])

            index_count = index_count + 1

            if(index_count > sim_dust_particles[size_bin]-1):
                index_count = 0
                size_bin = size_bin + 1


        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.scatter(r_p[:,0], r_p[:,1], r_p[:,2], s = 1)
        ax.set_title("Time = 0 s")
        ax.axes.set_xlim3d(left=1900, right=2100) 
        ax.axes.set_ylim3d(bottom=0, top=8) 
        ax.axes.set_zlim3d(bottom=0, top=1) 
        
        
        print("r_p at t") 
        print(r_p) 
        print("v_p at t") 
        print(v_p) 
        plt.savefig("monte_carlo_ejecta_"+str(j)+"_"+"0.png")
        plt.show()


        g_p = np.ones(int(sum(sim_dust_particles))) * 1.62
        

        # loop through particle move time
        dt_sub = 0.05
        for t in range(0,50):
            for p in range(0, int(sum(sim_dust_particles))):
                

                r_p[p,0] = r_p[p,0] + v_p[p,0] * dt_sub 
                r_p[p,1] = r_p[p,1] + v_p[p,1] * dt_sub 
                r_p[p,2] = r_p[p,2] + v_p[p,2] * dt_sub - 0.5 * g_p[p] * dt_sub**2

                v_p[p,0] = v_p[p,0]
                v_p[p,1] = v_p[p,1]
                v_p[p,2] = v_p[p,2] - g_p[p] * dt_sub


                if(r_p[p,2] <= 0):
                    v_p[p,:] = 0
                    g_p[p] = 0

            fig = plt.figure()
            ax = fig.add_subplot(projection='3d')
            ax.scatter(r_p[:,0], r_p[:,1], r_p[:,2], s = 1)
            ax.set_title("Time = "+str('{:.2f}'.format(dt_sub*(t+1)))+" s")
            ax.axes.set_xlim3d(left=1900, right=2100) 
            ax.axes.set_ylim3d(bottom=0, top=8) 
            ax.axes.set_zlim3d(bottom=0, top=1) 
            plt.savefig("monte_carlo_ejecta_"+str(j)+"_"+str(t+1)+".png")

            print("r_p at t") 
            print(r_p) 
            print("v_p at t") 
            print(v_p)  
            plt.close()

    
    #print(ejecta_mass_excavated)
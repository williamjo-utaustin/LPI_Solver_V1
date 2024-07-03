import numpy as np
import var_soil as soil


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

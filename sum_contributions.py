import numpy as np
import matplotlib.pyplot as plt


index_grain_sizes_list = np.zeros([4,2])

index_grain_sizes_list[0,:] = np.array([int(0), int(9)])
index_grain_sizes_list[1,:] = np.array([int(9), int(18)])
index_grain_sizes_list[2,:] = np.array([int(18),int(27)])
index_grain_sizes_list[3,:] = np.array([int(27),int(36)])

label_list = np.array(["1 - 10 microns", "10 - 90 microns", "100 - 900 microns", "1,000 - 10,000 microns"])
r_d = np.logspace(-1,4,50000)

print(label_list[0])


sum_energy_interp = np.zeros_like(r_d)
sum_mass_interp = np.zeros_like(r_d)


plt.figure(figsize=(5,5))


for i in range(0,4):

    print(i) 
    index_grain_sizes = index_grain_sizes_list[i,:]

    range_bounds_mid = np.genfromtxt("output/range_bounds_"+str(int(index_grain_sizes[0]))+"_"+str(int(index_grain_sizes[1]))+".csv", delimiter=",")
    mass_flux_upon_impact = np.genfromtxt("output/mass_flux_upon_impact_"+str(int(index_grain_sizes[0]))+"_"+str(int(index_grain_sizes[1]))+".csv", delimiter=',')
    momentum_flux_upon_impact = np.genfromtxt("output/momentum_flux_upon_impact_"+str(int(index_grain_sizes[0]))+"_"+str(int(index_grain_sizes[1]))+".csv", delimiter=',')
    energy_flux_upon_impact = np.genfromtxt("output/energy_flux_upon_impact_"+str(int(index_grain_sizes[0]))+"_"+str(int(index_grain_sizes[1]))+".csv", delimiter=',')

    total_mass_flux_upon_impact = np.zeros(np.size(range_bounds_mid))
    total_momentum_flux_upon_impact = np.zeros(np.size(range_bounds_mid))
    total_energy_flux_upon_impact = np.zeros(np.size(range_bounds_mid))

    for j in range(0,np.size(range_bounds_mid)):
        total_mass_flux_upon_impact[j] = np.sum(mass_flux_upon_impact[:,j])
        total_momentum_flux_upon_impact[j] = np.sum(momentum_flux_upon_impact[:,j])
        total_energy_flux_upon_impact[j] = np.sum(energy_flux_upon_impact[:,j])


    energy_interp = np.interp(r_d, range_bounds_mid, total_energy_flux_upon_impact, right = 0)
    sum_energy_interp = sum_energy_interp + energy_interp
    
    mass_interp = np.interp(r_d, range_bounds_mid, total_mass_flux_upon_impact, right = 0)
    sum_mass_interp = sum_mass_interp + mass_interp
    
    
    
    plt.semilogy(r_d, energy_interp, linestyle = 'dashed', linewidth = 2, label = label_list[i])
    #plt.semilogy(r_d, mass_interp, linestyle = 'dashed', linewidth = 2, label = label_list[i], zorder = (i+100) * 2)

plt.semilogy(r_d, sum_energy_interp, color = 'black', linewidth = 2, label = "Total Contribution")
#plt.semilogy(r_d, sum_mass_interp, color = 'black', linewidth = 2, label = "Total Contribution", zorder = 0)

plt.grid()
plt.legend(fontsize = 10)
plt.xlabel("Distance from Centerline (m)", fontsize = 16)
#plt.ylabel("Impact Mass per Unit Area $(kg/m^2)$", fontsize = 16)
plt.ylabel("Impact Energy per Unit Area $(J/m^2)$", fontsize = 16)
plt.xticks(fontsize = 18)
plt.yticks(fontsize = 18)
plt.locator_params(axis='x', nbins=5)

plt.xlim(0,250)
plt.ylim(1E0,1E4)
#plt.savefig("output/impact_mass_profile_combined_near.png", bbox_inches='tight',dpi=100)
plt.savefig("output/impact_energy_profile_combined_near.png", bbox_inches='tight',dpi=100)

#plt.xlim(0,8000)
#plt.ylim(1E0,1E4)
#plt.savefig("output/impact_mass_profile_combined_far.png", bbox_inches='tight',dpi=100)
#plt.savefig("output/impact_energy_profile_combined_far.png", bbox_inches='tight',dpi=100)

plt.show()
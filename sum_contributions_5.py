import numpy as np
import matplotlib.pyplot as plt

label_list = np.array(["1 - 10 $\mu m$", "10 - 90 $\mu m$", "100 - 900 $\mu m$", "1,000 - 10,000 $\mu m$"])
r_d = np.logspace(-1,5,50000)

index_grain_sizes_list = np.zeros([4,2])

index_grain_sizes_list[0,:] = np.array([int(0), int(9)])
index_grain_sizes_list[1,:] = np.array([int(9), int(18)])
index_grain_sizes_list[2,:] = np.array([int(18),int(27)])
index_grain_sizes_list[3,:] = np.array([int(27),int(36)])

sum_energy_interp = np.zeros_like(r_d)
sum_mass_interp = np.zeros_like(r_d)
sum_count_interp = np.zeros_like(r_d)


#file_folder = 'starship_raptor_nominal_50'
#identifier = '_scaled_scaling_9'

#file_folder = 'apollo_vel_slow_alt_mid'
#identifier = '_4_scaling_9'

file_folder = 'bluemoon_H30_lambda1p5_m1p5e04'
identifier = '_scaling_9'

#file_folder = 'starship_thruster_vel_high_50'
#file_folder = 'starship_thruster_nominal_100'
#file_folder = 'apollo_vel_slow_alt_high'
#file_folder = 'apollo_vel_slow_alt_mid'
#file_folder = 'apollo_vel_fast_alt_mid'

#ile_folder_array = np.array(['apollo_vel_slow_alt_mid', 'apollo_vel_mid_alt_mid', 'apollo_vel_fast_alt_mid'])


print(file_folder+identifier+"/"+file_folder+"_post_processing/range_bounds.csv")
range_bounds_mid = np.genfromtxt(file_folder+identifier+"/"+file_folder+"_post_processing/range_bounds.csv", delimiter=",")
mass_flux_upon_impact = np.genfromtxt(file_folder+identifier+"/"+file_folder+"_post_processing/"+file_folder+"_mass_flux_upon_impact.csv", delimiter=',')
energy_flux_upon_impact = np.genfromtxt(file_folder+identifier+"/"+file_folder+"_post_processing/"+file_folder+"_energy_flux_upon_impact.csv", delimiter=',')
count_flux_upon_impact = np.genfromtxt(file_folder+identifier+"/"+file_folder+"_post_processing/"+file_folder+"_count_flux_upon_impact.csv", delimiter=',')
        
count_energy_upon_impact_1 = np.genfromtxt(file_folder+identifier+"/"+file_folder+"_post_processing/"+file_folder+"_count_energy_threshold_impact_1.csv", delimiter=",")
count_energy_upon_impact_2 = np.genfromtxt(file_folder+identifier+"/"+file_folder+"_post_processing/"+file_folder+"_count_energy_threshold_impact_2.csv", delimiter=",")
count_energy_upon_impact_3 = np.genfromtxt(file_folder+identifier+"/"+file_folder+"_post_processing/"+file_folder+"_count_energy_threshold_impact_3.csv", delimiter=",")
count_energy_upon_impact_4 = np.genfromtxt(file_folder+identifier+"/"+file_folder+"_post_processing/"+file_folder+"_count_energy_threshold_impact_4.csv", delimiter=",")
count_energy_upon_impact_5 = np.genfromtxt(file_folder+identifier+"/"+file_folder+"_post_processing/"+file_folder+"_count_energy_threshold_impact_5.csv", delimiter=",")
count_energy_upon_impact_6 = np.genfromtxt(file_folder+identifier+"/"+file_folder+"_post_processing/"+file_folder+"_count_energy_threshold_impact_6.csv", delimiter=",")
count_energy_upon_impact_7 = np.genfromtxt(file_folder+identifier+"/"+file_folder+"_post_processing/"+file_folder+"_count_energy_threshold_impact_7.csv", delimiter=",")

plotting_variable = np.array(["count", "mass", "energy", "sample", "count_energy_1", "count_energy_2", "count_energy_3", "count_energy_4", "count_energy_5", "count_energy_6", "count_energy_7"])
linestyle_array = np.array(['solid', 'dotted', 'dashed', 'dashdot'])

# loop through count, mass, energy, and sample
for n in range(0,11):
    plt.figure(figsize=(5,6))

    to_plot = plotting_variable[n]
    print("Printing ", to_plot)

    total_mass_flux_upon_impact = np.zeros(np.size(range_bounds_mid))
    total_energy_flux_upon_impact = np.zeros(np.size(range_bounds_mid))
    total_count_flux_upon_impact = np.zeros(np.size(range_bounds_mid))

    total_count_energy_upon_impact_1 = np.zeros(np.size(range_bounds_mid))
    total_count_energy_upon_impact_2 = np.zeros(np.size(range_bounds_mid))
    total_count_energy_upon_impact_3 = np.zeros(np.size(range_bounds_mid))
    total_count_energy_upon_impact_4 = np.zeros(np.size(range_bounds_mid))
    total_count_energy_upon_impact_5 = np.zeros(np.size(range_bounds_mid))
    total_count_energy_upon_impact_6 = np.zeros(np.size(range_bounds_mid))
    total_count_energy_upon_impact_7 = np.zeros(np.size(range_bounds_mid))

    component_mass_flux_upon_impact = np.zeros(np.size(range_bounds_mid))
    component_energy_flux_upon_impact = np.zeros(np.size(range_bounds_mid))
    component_count_flux_upon_impact = np.zeros(np.size(range_bounds_mid))
    
    component_count_energy_upon_impact_1 = np.zeros(np.size(range_bounds_mid))
    component_count_energy_upon_impact_2 = np.zeros(np.size(range_bounds_mid))
    component_count_energy_upon_impact_3 = np.zeros(np.size(range_bounds_mid))
    component_count_energy_upon_impact_4 = np.zeros(np.size(range_bounds_mid))
    component_count_energy_upon_impact_5 = np.zeros(np.size(range_bounds_mid))
    component_count_energy_upon_impact_6 = np.zeros(np.size(range_bounds_mid))
    component_count_energy_upon_impact_7 = np.zeros(np.size(range_bounds_mid))

    for j in range(0,np.size(range_bounds_mid)):

        total_mass_flux_upon_impact[j] = np.sum(mass_flux_upon_impact[:,j])
        total_energy_flux_upon_impact[j] = np.sum(energy_flux_upon_impact[:,j])
        total_count_flux_upon_impact[j] = np.sum(count_flux_upon_impact[:,j])
    
        total_count_energy_upon_impact_1[j] = np.sum(count_energy_upon_impact_1[:,j])
        total_count_energy_upon_impact_2[j] = np.sum(count_energy_upon_impact_2[:,j])
        total_count_energy_upon_impact_3[j] = np.sum(count_energy_upon_impact_3[:,j])
        total_count_energy_upon_impact_4[j] = np.sum(count_energy_upon_impact_4[:,j])
        total_count_energy_upon_impact_5[j] = np.sum(count_energy_upon_impact_5[:,j])
        total_count_energy_upon_impact_6[j] = np.sum(count_energy_upon_impact_6[:,j])
        total_count_energy_upon_impact_7[j] = np.sum(count_energy_upon_impact_7[:,j])


    for i in range(0,4):
        for j in range(0,np.size(range_bounds_mid)):
            component_mass_flux_upon_impact[j] = np.sum(mass_flux_upon_impact[int(index_grain_sizes_list[i,0]):int(index_grain_sizes_list[i,1]),j])
            component_energy_flux_upon_impact[j] = np.sum(energy_flux_upon_impact[int(index_grain_sizes_list[i,0]):int(index_grain_sizes_list[i,1]),j])
            component_count_flux_upon_impact[j] = np.sum(count_flux_upon_impact[int(index_grain_sizes_list[i,0]):int(index_grain_sizes_list[i,1]),j])

            component_count_energy_upon_impact_1[j] = np.sum(count_energy_upon_impact_1[int(index_grain_sizes_list[i,0]):int(index_grain_sizes_list[i,1]),j])
            component_count_energy_upon_impact_2[j] = np.sum(count_energy_upon_impact_2[int(index_grain_sizes_list[i,0]):int(index_grain_sizes_list[i,1]),j])
            component_count_energy_upon_impact_3[j] = np.sum(count_energy_upon_impact_3[int(index_grain_sizes_list[i,0]):int(index_grain_sizes_list[i,1]),j])
            component_count_energy_upon_impact_4[j] = np.sum(count_energy_upon_impact_4[int(index_grain_sizes_list[i,0]):int(index_grain_sizes_list[i,1]),j])
            component_count_energy_upon_impact_5[j] = np.sum(count_energy_upon_impact_5[int(index_grain_sizes_list[i,0]):int(index_grain_sizes_list[i,1]),j])
            component_count_energy_upon_impact_6[j] = np.sum(count_energy_upon_impact_6[int(index_grain_sizes_list[i,0]):int(index_grain_sizes_list[i,1]),j])
            component_count_energy_upon_impact_7[j] = np.sum(count_energy_upon_impact_7[int(index_grain_sizes_list[i,0]):int(index_grain_sizes_list[i,1]),j])

        #plt.semilogy(range_bounds_mid, component_count_flux_upon_impact * 0.0025, linestyle = 'dashed', label = label_list[i])
        #plt.semilogy(range_bounds_mid, component_energy_flux_upon_impact, linestyle = 'dashed', label = label_list[i])

        if(to_plot == 'mass'):
            plt.loglog(range_bounds_mid, component_mass_flux_upon_impact, linestyle = linestyle_array[i], label = label_list[i], zorder = i+1, linewidth = 3)
        elif(to_plot == 'energy'):
            plt.loglog(range_bounds_mid, component_energy_flux_upon_impact, linestyle = linestyle_array[i], label = label_list[i], zorder = i+1, linewidth = 3)
        elif(to_plot == 'count'):
            plt.loglog(range_bounds_mid, component_count_flux_upon_impact, linestyle = linestyle_array[i], label = label_list[i], zorder = i+1, linewidth = 3)
        elif(to_plot == 'sample'):
            plt.loglog(range_bounds_mid, component_count_flux_upon_impact * 0.0025, linestyle = linestyle_array[i], label = label_list[i], zorder = i+1, linewidth = 3)
        elif(to_plot == 'count_energy_1'):
            plt.loglog(range_bounds_mid, component_count_energy_upon_impact_1, linestyle = linestyle_array[i], label = label_list[i], zorder = i+1, linewidth = 3)
        elif(to_plot == 'count_energy_2'):
            plt.loglog(range_bounds_mid, component_count_energy_upon_impact_2, linestyle = linestyle_array[i], label = label_list[i], zorder = i+1, linewidth = 3)
        elif(to_plot == 'count_energy_3'):
            plt.loglog(range_bounds_mid, component_count_energy_upon_impact_3, linestyle = linestyle_array[i], label = label_list[i], zorder = i+1, linewidth = 3)
        elif(to_plot == 'count_energy_4'):
            plt.loglog(range_bounds_mid, component_count_energy_upon_impact_4, linestyle = linestyle_array[i], label = label_list[i], zorder = i+1, linewidth = 3)
        elif(to_plot == 'count_energy_5'):
            plt.loglog(range_bounds_mid, component_count_energy_upon_impact_5, linestyle = linestyle_array[i], label = label_list[i], zorder = i+1, linewidth = 3)
        elif(to_plot == 'count_energy_6'):
            plt.loglog(range_bounds_mid, component_count_energy_upon_impact_6, linestyle = linestyle_array[i], label = label_list[i], zorder = i+1, linewidth = 3)
        elif(to_plot == 'count_energy_7'):
            plt.loglog(range_bounds_mid, component_count_energy_upon_impact_7, linestyle = linestyle_array[i], label = label_list[i], zorder = i+1, linewidth = 3)
        else:
            exit("Unknown Variable Plotted")

    #plt.semilogy(range_bounds_mid, total_count_flux_upon_impact * 0.0025, color = 'black', linewidth = 2, label = 'Total Contribution')
    if(to_plot == 'mass'):
        plt.ylabel("Impact Mass $(kg/m^2)$", fontsize = 16)
        plt.loglog(range_bounds_mid, total_mass_flux_upon_impact, color = 'black', linewidth = 4, label = 'Total Contribution', zorder = 0)
    elif(to_plot == 'energy'):
        plt.ylabel("Impact Energy $(J/m^2)$", fontsize = 16)
        plt.loglog(range_bounds_mid, total_energy_flux_upon_impact, color = 'black', linewidth = 4, label = 'Total Contribution', zorder = 0)
    elif(to_plot == 'count'):
        plt.ylabel("Impact Count $(\#/m^2)$", fontsize = 16)
        plt.loglog(range_bounds_mid, total_count_flux_upon_impact, color = 'black', linewidth = 4, label = 'Total Contribution', zorder = 0)
    elif(to_plot == 'sample'):
        plt.ylabel("Impact Count $(\#/25 cm^2)$", fontsize = 16)
        plt.loglog(range_bounds_mid, total_count_flux_upon_impact * 0.0025, color = 'black', linewidth = 4, label = 'Total Contribution', zorder = 0)
    elif(to_plot == 'count_energy_1'):
        plt.ylabel("Impact Count $(\#/m^2)$", fontsize = 16)
        plt.loglog(range_bounds_mid, total_count_energy_upon_impact_1, color = 'black', linewidth = 4, label = 'Total Contribution', zorder = 0)
    elif(to_plot == 'count_energy_2'):
        plt.ylabel("Impact Count $(\#/m^2)$", fontsize = 16)
        plt.loglog(range_bounds_mid, total_count_energy_upon_impact_2, color = 'black', linewidth = 4, label = 'Total Contribution', zorder = 0)
    elif(to_plot == 'count_energy_3'):
        plt.ylabel("Impact Count $(\#/m^2)$", fontsize = 16)
        plt.loglog(range_bounds_mid, total_count_energy_upon_impact_3, color = 'black', linewidth = 4, label = 'Total Contribution', zorder = 0)
    elif(to_plot == 'count_energy_4'):
        plt.ylabel("Impact Count $(\#/m^2)$", fontsize = 16)
        plt.loglog(range_bounds_mid, total_count_energy_upon_impact_4, color = 'black', linewidth = 4, label = 'Total Contribution', zorder = 0)
    elif(to_plot == 'count_energy_5'):
        plt.ylabel("Impact Count $(\#/m^2)$", fontsize = 16)
        plt.loglog(range_bounds_mid, total_count_energy_upon_impact_5, color = 'black', linewidth = 4, label = 'Total Contribution', zorder = 0)
    elif(to_plot == 'count_energy_6'):
        plt.ylabel("Impact Count $(\#/m^2)$", fontsize = 16)
        plt.loglog(range_bounds_mid, total_count_energy_upon_impact_6, color = 'black', linewidth = 4, label = 'Total Contribution', zorder = 0)
    elif(to_plot == 'count_energy_7'):
        plt.ylabel("Impact Count $(\#/m^2)$", fontsize = 16)
        plt.loglog(range_bounds_mid, total_count_energy_upon_impact_7, color = 'black', linewidth = 4, label = 'Total Contribution', zorder = 0)
    else:
        exit("Unknown Variable Plotted")


    legend = plt.legend(fontsize = 12, loc = 'upper right')
    legend.set_title("Particle Diameters", prop={'size':12}) 
    #plt.locator_params(axis='x', nbins=5)
    plt.grid()
    #plt.xlim(1E0, 1E8)
    plt.xlim(1E0, 1E10)
    #plt.ylim(1E0,1E13)
    plt.xticks(fontsize = 18)
    plt.yticks(fontsize = 18)
    plt.xlabel("Distance from Plume" + "\n" + "Centerline (m)", fontsize = 16)
    #plt.ylabel("Impact Count per Sample Size $(\#/25cm^2)$", fontsize = 14)
    #plt.savefig("output_starship_vel_mid/impact_count_profile_combined_near.png", bbox_inches='tight',dpi=100)

    if(to_plot == 'mass'):
        #plt.ylim(1E-14,1E2)
        plt.ylim(1E-10,1E0)
        plt.savefig(file_folder+identifier+"/"+file_folder+"_post_processing/mass.png", bbox_inches='tight',dpi=100)
    elif(to_plot == 'energy'):
        #plt.ylim(1E-8,1E4)
        plt.ylim(1E-8,1E2)
        plt.savefig(file_folder+identifier+"/"+file_folder+"_post_processing/energy.png", bbox_inches='tight',dpi=100)
    elif(to_plot == 'count'):
        #plt.ylim(1E0,1E14)
        plt.ylim(1E-2,1E12)
        plt.savefig(file_folder+identifier+"/"+file_folder+"_post_processing/count.png", bbox_inches='tight',dpi=100)
    elif(to_plot == 'sample'):
        #plt.ylim(1E-3,1E11)
        plt.ylim(1E-2,1E10)
        plt.savefig(file_folder+identifier+"/"+file_folder+"_post_processing/sample_count.png", bbox_inches='tight',dpi=100)
    elif(to_plot == 'count_energy_1'):
        #plt.ylim(1E0,1E14)
        plt.ylim(1E-2,1E12)
        plt.savefig(file_folder+identifier+"/"+file_folder+"_post_processing/count_energy_1.png", bbox_inches='tight',dpi=100)
    elif(to_plot == 'count_energy_2'):
        #plt.ylim(1E0,1E14)
        plt.ylim(1E-2,1E12)
        plt.savefig(file_folder+identifier+"/"+file_folder+"_post_processing/count_energy_2.png", bbox_inches='tight',dpi=100)
    elif(to_plot == 'count_energy_3'):
        #plt.ylim(1E0,1E14)
        plt.ylim(1E-2,1E12)
        plt.savefig(file_folder+identifier+"/"+file_folder+"_post_processing/count_energy_3.png", bbox_inches='tight',dpi=100)
    elif(to_plot == 'count_energy_4'):
        #plt.ylim(1E0,1E14)
        plt.ylim(1E-2,1E12)
        plt.savefig(file_folder+identifier+"/"+file_folder+"_post_processing/count_energy_4.png", bbox_inches='tight',dpi=100)
    elif(to_plot == 'count_energy_5'):
        #plt.ylim(1E0,1E14)
        plt.ylim(1E-2,1E12)
        plt.savefig(file_folder+identifier+"/"+file_folder+"_post_processing/count_energy_5.png", bbox_inches='tight',dpi=100)
    elif(to_plot == 'count_energy_6'):
        #plt.ylim(1E0,1E14)
        plt.ylim(1E-2,1E12)
        plt.savefig(file_folder+identifier+"/"+file_folder+"_post_processing/count_energy_6.png", bbox_inches='tight',dpi=100)
    elif(to_plot == 'count_energy_7'):
        #plt.ylim(1E0,1E14)
        plt.ylim(1E-2,1E12)
        plt.savefig(file_folder+identifier+"/"+file_folder+"_post_processing/count_energy_7.png", bbox_inches='tight',dpi=100)
    
    else:
        exit("Unknown Variable Plotted")

    plt.close()

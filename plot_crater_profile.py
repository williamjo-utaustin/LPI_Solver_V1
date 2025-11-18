import numpy as np
import matplotlib.pyplot as plt

# this program plots the final crater size


#folder = 'apollo_vel_fast_alt_mid' #287
#folder = 'apollo_vel_mid_alt_mid' #526
#folder = 'apollo_vel_slow_alt_mid' #3150

#folder_array = np.array(['apollo_vel_fast_alt_mid', 'apollo_vel_mid_alt_mid', 'apollo_vel_slow_alt_mid'])
#end_index_array = np.array([287, 526, 3150])
#label_array = np.array(['-11 m/s', '-6 m/s', '-1 m/s'])
#color_array = np.array(['tab:blue', 'tab:orange', 'tab:green'])
#linestyle_array = np.array(['solid', 'dotted', 'dashed'])


#folder_array = np.array(['apollo_vel_slow_alt_high', 'apollo_vel_slow_alt_mid', 'apollo_vel_slow_alt_low'])
#end_index_array = np.array([4401, 3150, 1401])
#label_array = np.array(['45 m', '32.5 m', '15 m'])
#color_array = np.array(['tab:purple', 'tab:green', 'tab:pink'])
#linestyle_array = np.array(['dashdot', 'dashed', 'solid'])

#folder_array = np.array(['apollo_vel_slow_alt_high', 'apollo_vel_slow_alt_mid', 'apollo_vel_slow_alt_low'])
#end_index_array = np.array([4401, 3150, 1401])
#label_array = np.array(['45 m', '32.5 m', '15 m'])
#index_x = np.array([75, 42, 51])

folder_array = np.array(['starship_thruster_vel_high_100', 'starship_thruster_nominal_100'])
#folder_array = np.array(['starship_thruster_vel_high_100', 'starship_raptor_vel_high_100'])
#folder_array = np.array(['starship_thruster_vel_high_50', 'starship_raptor_vel_high_50'])
end_index_array = np.array([2273, 8801])
label_array = np.array(['Upper Thrusters (UT)', 'Lower Thrusters (RV)'])
color_array = np.array(['tab:blue', 'tab:orange'])
linestyle_array = np.array(['solid', 'dotted'])
#index_x = np.array([75, 42, 51])


plt.figure(figsize=(12,4))

for f in range(0,2):
    folder = folder_array[f]
    end_index = str(end_index_array[f])

    for i in range(0, int(end_index)):

        timestep = str(i)
        data_file = folder +"/"+folder+"_ejecta_properties_"+timestep+".csv"
        data = np.genfromtxt(data_file, delimiter = ' ')
        
        if(i == 0):

            # obtain the final array size and stuff into array
            final_array_point = int(np.size(data[:,0]))-1
            array_size = int(data[final_array_point, 0])+1
            save_array_init_x = np.zeros(array_size)
            save_array_init_y = np.zeros(array_size)

            for j in range(0, int(np.size(data[:,0]))):
                save_array_init_x[int(data[j,0])] = data[j,1] 
                save_array_init_y[int(data[j,0])] = data[j,4]

        else:
            save_array = data[:,4]
            for j in range(0, np.size(save_array)):
                save_array_init_x[int(data[j,0])] = data[j,1] 
                save_array_init_y[int(data[j,0])] = data[j,4]

    #save_array_init_x = np.append(save_array_init_x, index_x[f])
    #save_array_init_y = np.append(save_array_init_y, 0)


    plt.plot(save_array_init_x, save_array_init_y, linewidth = 5, label = label_array[f], color = color_array[f], linestyle = linestyle_array[f])

plt.xlim(0, 750)
plt.ylim(0.6, 0)
#plt.xlim(0, 50)
#plt.ylim(0.6, 0)
plt.grid()
plt.xticks(fontsize = 18)
plt.yticks(fontsize = 18)
plt.ylabel("Excavation Depth (m)", fontsize = 18)
plt.xlabel("Distance from Plume Centerline (m)", fontsize = 18)
legend = plt.legend(fontsize = 14)
legend.set_title('Thruster Configuration', prop={'size':14})
plt.savefig("final_crater_profile_together_1.png", bbox_inches='tight',dpi=100)
print(save_array)

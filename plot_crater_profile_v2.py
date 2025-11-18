import numpy as np
import matplotlib.pyplot as plt

# this program plots the final crater size

plot_on = False


#num_curves_per_plot = 1
#folder_array = np.array(['starship_raptor_nominal_50'])
#identifier_array = np.array(['_scaled_dec2024_scaling_9'])
#end_index_array = np.array([867])
#label_array = np.array(['HLS Lower Thrusters - 50 Tons'])
#color_array = np.array(['tab:blue'])
#solid_array = np.array(['solid'])
#linestyle_array = np.array(['dashed'])
#end_range = 100
#time = 0
#altitude = 31.5

#starship_raptor_nominal_100_scaled_dec2024_scaling_9/

num_curves_per_plot = 1
folder_array = np.array(['apollo_vel_slow_alt_mid'])
identifier_array = np.array(['_4_scaling_9'])
end_index_array = np.array([315])
label_array = np.array(['Apollo LM'])
color_array = np.array(['tab:blue'])
solid_array = np.array(['solid'])
end_range = 100
time = 0
altitude = 31.5

#num_curves_per_plot = 2
#folder_array = np.array(['starship_thruster_nominal_100', 'starship_raptor_nominal_100'])
#identifier_array = np.array(['_scaled_scaling_9', '_scaled_scaling_9'])
#end_index_array = np.array([886, 886])
#label_array = np.array(['HLS Upper Thrusters - 50 Tons', 'HLS Lower Thrusters - 50 Tons'])
#color_array = np.array(['tab:orange', 'tab:green'])
#linestyle_array = np.array(['dashed', 'dotted'])
#end_range = 100
#time = 0
#altitude = 250

for t in range(311,int(end_index_array[0])+1,1):
   
    if(plot_on == True):
        plt.figure(figsize=(12,4))
    
    time = time + 0.1
    
    # at each timestep, we open both folders or many folders
    for f in range(0,num_curves_per_plot):
        
        folder = folder_array[f]
        identifier = identifier_array[f]
        end_index = str(end_index_array[f])

        # collect data from previous timesteps, as our data is 'concatenated' to save space
        for i in range(0, t):

            
            timestep = str(i)
            data_file = folder+identifier +"/"+folder+"_ejecta_properties_"+timestep+".csv"
            data = np.genfromtxt(data_file, delimiter = ' ')
            
            # create the array size from the first timestep (will all be zero)
            if(i == 0):

                # obtain the final array size and stuff into array
                final_array_point = int(np.size(data[:,0]))-1
                array_size = int(data[final_array_point, 0])+1
                save_array_init_x = np.zeros([array_size,num_curves_per_plot])
                save_array_init_y = np.zeros([array_size,num_curves_per_plot])

                for j in range(0, int(np.size(data[:,0]))):
                    save_array_init_x[int(data[j,0]), f] = data[j,1] 
                    save_array_init_y[int(data[j,0]), f] = data[j,4]

            # alter the array
            else:
                save_array = data[:,4]
                for j in range(0, np.size(save_array)):
                    save_array_init_x[int(data[j,0]), f] = data[j,1] 
                    save_array_init_y[int(data[j,0]), f] = data[j,4]

        #save_array_init_x = np.append(save_array_init_x, index_x[f])
        #save_array_init_y = np.append(save_array_init_y, 0)

        # once previous timesteps have been collected, plot the profile
        if(plot_on == True):
            plt.plot(save_array_init_x[:,f], save_array_init_y[:,f]*100, linewidth = 5, label = label_array[f], color = color_array[f], linestyle = linestyle_array[f])
       

        if(t == int(end_index_array[0]) - 1 or np.mod(t,100) == 0):
            print("Plotting Profile ", t)
            
        x = save_array_init_x[:, f]
        y = save_array_init_y[:, f] * 100
        
        idx_max = np.argmax(y)         # index of max y
        x_at_max_y = x[idx_max]        # corresponding x value
        y_max = y[idx_max]             # (optional) max y value
        
        print(t,",",f"Max y = {y_max:.3f} occurs at x = {x_at_max_y:.3f}")





    
    if(altitude >= 200):
        v_desc = -11
    elif(altitude < 200 and altitude >= 150):
        v_desc = -8.5
    elif(altitude < 150 and altitude >= 100):
        v_desc = -6
    elif(altitude < 100 and altitude >= 50):
        v_desc = -2.5
    elif(altitude < 50 and altitude >= 0):
        v_desc = -1
    else:
        v_desc = -1

    altitude = altitude + 0.1 * v_desc

    






    if(plot_on == True):

        plt.title("Time = "+str("{:.2f}".format((time)))+" s" + ", Lander Altitude = "+str("{:.2f}".format(altitude))+" m"+ ", Descent Velocity = "+str("{:.2f}".format(v_desc)) + " m/s", fontsize = 18)
        plt.xlim(0, end_range)
        plt.ylim(6, 0)
        #plt.xlim(0, 50)
        #plt.ylim(0.6, 0)
        plt.grid()
        plt.xticks(fontsize = 18)
        plt.yticks(fontsize = 18)
        plt.ylabel("Excavation Depth (cm)", fontsize = 18)
        plt.xlabel("Distance from Plume Centerline (m)", fontsize = 18)
        legend = plt.legend(fontsize = 14, loc = 'lower right')
        legend.set_title('Thruster Configuration', prop={'size':14})
        plt.savefig(folder+identifier +"/crater_profile"+"_"+str(t)+".png", bbox_inches='tight',dpi=100)
        print(folder+identifier +"/crater_profile"+"_"+str(t)+".png")
        
        
        if (t == int(end_index_array[0])-1):
            plt.title("Time = "+str("{:.2f}".format((t+1)*0.1))+" s" + ", Lander Altitude = "+str("{:.2f}".format(0))+" m" + ", Descent Velocity = "+str("{:.2f}".format(0)) + " m/s", fontsize = 18)
            plt.savefig(folder+identifier +"/crater_profile"+"_"+str(t+1)+".png", bbox_inches='tight',dpi=100)
             
        plt.close()
        #print(save_array)

        # plot the initial condition
        if (t == 1):
            
            plt.figure(figsize=(12,4))
            
            for f in range(0,num_curves_per_plot):
                plt.plot(save_array_init_x[:,f], save_array_init_y[:,f]*0, linewidth = 5, label = label_array[f], color = color_array[f], linestyle = linestyle_array[f])
            
            plt.title("Time = "+str("{:.2f}".format((0)*0.1))+" s" + ", Lander Altitude = "+str("{:.2f}".format(250))+" m" + ", Descent Velocity = "+str("{:.2f}".format(-11)) + " m/s", fontsize = 18)
            plt.xlim(0, end_range)
            plt.ylim(5, 0)
            plt.grid()
            plt.xticks(fontsize = 18)
            plt.yticks(fontsize = 18)
            plt.ylabel("Excavation Depth (cm)", fontsize = 18)
            plt.xlabel("Distance from Plume Centerline (m)", fontsize = 18)
            legend = plt.legend(fontsize = 14, loc = 'lower right')
            legend.set_title('Thruster Configuration', prop={'size':14})
            plt.savefig(folder+identifier +"/crater_profile"+"_"+str(0)+".png", bbox_inches='tight',dpi=100)
            plt.close()

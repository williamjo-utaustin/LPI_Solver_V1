import numpy as np
import matplotlib.pyplot as plt
import os
# this program plots the final crater size


dt = 0.01
starting_altitude = 30
v_init = -19.718
a_upward = 6.48
t_td = 15

n_timesteps = int(t_td/dt)

plot_on = True

num_curves_per_plot = 1
#folder_array = np.array(['starship_H60_lambda10_m1e05'])
folder_array = np.array(['starship_H60_lambda10_m5e04'])
#folder_array = np.array(['bluemoon_H60_lambda4_m1p5e04'])
identifier_array = np.array(['_scaling_9'])
end_index_array = np.array([n_timesteps])
label_array = np.array(['Bluemoon'])
#label_array = np.array(['Blue Origin Test'])
color_array = np.array(['tab:blue'])
solid_array = np.array(['solid'])
linestyle_array = np.array(['solid'])
end_range = 25
time = 0


for t in range(280,int(end_index_array[0]), 1):
   
    if(plot_on == True):
        plt.figure(figsize=(12,4))
    
    time = time + dt
    
    # at each timestep, we open both folders or many folders
    for f in range(0,num_curves_per_plot):
        
        folder = folder_array[f]
        identifier = identifier_array[f]
        end_index = str(end_index_array[f])
        output_dir = folder + identifier
        os.makedirs(output_dir, exist_ok=True)

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
            
            #for j in range(0, array_size):
            #    print("Curve ", f, save_array_init_x[j,f], save_array_init_y[j,f])
        

    if(plot_on == True):

        #plt.title("Time = "+str("{:.2f}".format((time)))+" s" + ", Lander Altitude = "+str("{:.2f}".format(altitude))+" m"+ ", Descent Velocity = "+str("{:.2f}".format(v_desc)) + " m/s", fontsize = 18)
        plt.xlim(0, end_range)
        plt.ylim(10, 0)
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

        x = save_array_init_x[:, 0]
        y = save_array_init_y[:, 0] * 100
        
        # Save to CSV
        csv_filename = f"{output_dir}/{output_dir}_crater_profile_{t}.csv"
        np.savetxt(csv_filename, np.column_stack((x, y)), delimiter=",", header="x,y", comments='')

        

        # plot the initial condition
        if (t == 1):
            
            plt.figure(figsize=(12,4))
            
            for f in range(0,num_curves_per_plot):
                plt.plot(save_array_init_x[:,f], save_array_init_y[:,f]*0, linewidth = 5, label = label_array[f], color = color_array[f], linestyle = linestyle_array[f])
            
            #plt.title("Time = "+str("{:.2f}".format((0)*0.1))+" s" + ", Lander Altitude = "+str("{:.2f}".format(250))+" m" + ", Descent Velocity = "+str("{:.2f}".format(-11)) + " m/s", fontsize = 18)
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

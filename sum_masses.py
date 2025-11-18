import numpy as np
import matplotlib.pyplot as plt

#file_foler = file_folder = 'starship_thruster_vel_high_100'
#file_foler = file_folder = 'starship_raptor_vel_high_100'


#file_folder = 'starship_H30_lambda1p5_m5e04'
#file_folder = 'starship_H30_lambda2_m5e04'
#file_folder = 'starship_H30_lambda4_m5e04'
#file_folder = 'starship_H30_lambda6_m5e04'
#file_folder = 'starship_H30_lambda8_m5e04'
#file_folder = 'starship_H30_lambda10_m5e04'
#
#file_folder = 'starship_H30_lambda1p5_m1e05'
#file_folder = 'starship_H30_lambda2_m1e05'
#file_folder = 'starship_H30_lambda4_m1e05'
#file_folder = 'starship_H30_lambda6_m1e05'
#file_folder = 'starship_H30_lambda8_m1e05'
#file_folder = 'starship_H30_lambda10_m1e05'
#
#file_folder = 'starship_H60_lambda1p5_m5e04'
#file_folder = 'starship_H60_lambda2_m5e04'
#file_folder = 'starship_H60_lambda4_m5e04'
#file_folder = 'starship_H60_lambda6_m5e04'
#file_folder = 'starship_H60_lambda8_m5e04'
#file_folder = 'starship_H60_lambda10_m5e04'
#
#file_folder = 'starship_H60_lambda1p5_m1e05'
#file_folder = 'starship_H60_lambda2_m1e05'
#file_folder = 'starship_H60_lambda4_m1e05'
#file_folder = 'starship_H60_lambda6_m1e05'
#file_folder = 'starship_H60_lambda8_m1e05'
#file_folder = 'starship_H60_lambda10_m1e05'




file_folder = 'bluemoon_H45_lambda2p5_m1p5e04'
#file_folder = 'bluemoon_H30_lambda4_m1p5e04'
#file_folder = 'apollo_vel_slow_alt_mid'
#file_folder = "starship_raptor_vel_high_50"
#file_folder = "starship_raptor_vel_high_100"
#file_folder = "starship_thruster_vel_high_50"
#file_folder = "starship_thruster_nominal_100"
#file_folder = "starship_thruster_vel_high_50"
#file_folder = "starship_raptor_nominal_100"

print(file_folder)
sum = 0
#plt.figure(figsize = (6,4))
#for t in range(0,867):
for t in range(0,10000):
#for t in range(0,287):

    # for new starship scaling
    #ejecta_properties = np.genfromtxt(file_folder+"_scaled_scaling_9/"+file_folder+"_ejecta_properties_"+str(t)+".csv", delimiter = " ")
    
    # for new apollo scaling
    #ejecta_properties = np.genfromtxt(file_folder+"_scaled_scaling_9/"+file_folder+"_ejecta_properties_"+str(t)+".csv", delimiter = " ")
    #ejecta_properties = np.genfromtxt(file_folder+"_scaled_dec2024_scaling_9/"+file_folder+"_ejecta_properties_"+str(t)+".csv", delimiter = " ")
    ejecta_properties = np.genfromtxt(file_folder+"_scaling_9/"+file_folder+"_ejecta_properties_"+str(t)+".csv", delimiter = " ")
    #ejecta_properties = np.genfromtxt(file_folder+"_4_scaling_9/"+file_folder+"_ejecta_properties_"+str(t)+".csv", delimiter = " ")
    #ejecta_properties = np.genfromtxt(file_folder+"/"+file_folder+"_ejecta_properties_"+str(t)+".csv", delimiter = " ")
    sum = sum +np.sum(ejecta_properties[:,3])
    print(t,",",sum)    
    #if(np.mod(t,5) == 0):
    #    plt.scatter(t*0.1, sum, color = 'tab:blue')
    
#plt.xlim(0,31.5)
#plt.ylim(0,3000)
#plt.xlabel("Time (s)", fontsize = 18)
#plt.ylabel("Mass Eroded (kg)", fontsize = 18)
#plt.xticks(fontsize = 18)
#plt.yticks(fontsize = 18)
#plt.grid()
#plt.savefig(file_folder+"_4_scaling_9/"+"mass_eroded_"+str(t)+".png", bbox_inches = 'tight', dpi = 100)
#plt.close()

print(sum/1000, "tons") 

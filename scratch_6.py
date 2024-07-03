import numpy as np
import matplotlib.pyplot as plt

#data = np.genfromtxt('output/ejecta_properties_0.csv', delimiter = ' ')
#
#for i in range(0,np.size(data[:,3])):
#    print(data[i,1], data[i,3])
#
#
#
#plt.figure(figsize=(5,5))
#plt.plot(data[:,1], data[:,3],        marker = 'o', label = "Total Contribution", linewidth = 2, color = 'black')
#plt.plot(data[:,1], data[:,3]*0.016,  marker = 'o', label = "Just 1 $\mu m$", linewidth = 2, linestyle = 'dashed')
#plt.plot(data[:,1], data[:,3]*0.1293, marker = 'o', label = "From 1 to 9 $\mu m$", linewidth = 2, linestyle = 'dashed')
#plt.plot(data[:,1], data[:,3]*0.4401, marker = 'o', label = "From 10 to 90 $\mu m$", linewidth = 2, linestyle = 'dashed')
#plt.plot(data[:,1], data[:,3]*0.3505, marker = 'o', label = "From 100 to 1000 $\mu m$", linewidth = 2, linestyle = 'dashed')
#plt.plot(data[:,1], data[:,3]*0.0801, marker = 'o', label = "From 1000 to 9000 $\mu m$", linewidth = 2, linestyle = 'dashed')
#
#plt.grid()
#plt.title("Mass Eroded between 0 - 0.01s", fontsize = 18)
#plt.legend(fontsize = 8)
#plt.xlabel("Ring Location from Centerline (m)", fontsize = 16)
#plt.ylabel("Mass Lofted by Ring Bin (kg)", fontsize = 16)
#plt.xticks(fontsize = 18)
#plt.yticks(fontsize = 18)
#plt.savefig("output/plot_mass_contributions.png", bbox_inches='tight',dpi=100)
#
#plt.show()

total_sum = 0
for i in range(0,301):
    data = np.genfromtxt('output/ejecta_properties_'+str(i)+'.csv', delimiter = ' ')
    total_sum = total_sum +np.sum(data[:,3]) 
    print(i, np.sum(data[:,3]))

print(total_sum)
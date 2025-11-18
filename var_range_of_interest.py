import numpy as np
import spacecraft_options as spacecraft

min_altitude = 1000
max_altitude = 5000
n_points_altitude = 9


if('apollo_' in spacecraft.lander_type):
    min_centerline = 0
    max_centerline = 50
    n_points_centerline = 50
    #max_centerline = 100
    #n_points_centerline = 25
    #n_points_centerline = 200

elif("starship_" in spacecraft.lander_type):
    min_centerline = 0
    max_centerline = 800
    n_points_centerline = 200 


elif("bluemoon_" in spacecraft.lander_type):
    min_centerline = 0
    max_centerline = 50
    n_points_centerline = 140 
    
else:
    print("Need to specify a valid lander type!")
    exit()



# option to set a 
h_nozzle_array = np.linspace(min_altitude, max_altitude, n_points_altitude)
#r_centerlines_array = np.linspace(min_centerline,max_centerline,n_points_centerline)
r_centerlines_array = np.geomspace(1E-1, max_centerline, num = n_points_centerline)

# 1,2,3,4,5,6,7,8,9,0, 2,3,4,5,6,7,8,9,0, 2,3,4,5,6,7,8,9,0, 2,3,4,5,6,7,8,9,0, 2,3,4,5,6,7,8,9,0

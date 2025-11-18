import numpy as np
import matplotlib.pyplot as plt

import var_constants as constants
import var_timestepping as timestep
import var_soil as soil
import sys_output as out
import var_range_of_interest as bounds
import var_impinged_gas as imp
import var_nozzle as nozzle
import spacecraft_options as spacecraft

import var_scale_erosion as scale

#import var_solve_ejection as sol_ej
if (spacecraft.lander_type == 'bluemoon_H30_lambda1p01_m1p5e04'):
    output_folder = 'bluemoon_H30_lambda1p01_m1p5e04'
elif (spacecraft.lander_type == 'bluemoon_H30_lambda1p01_m3e04'):
    output_folder = 'bluemoon_H30_lambda1p01_m3e04'
elif (spacecraft.lander_type == 'bluemoon_H30_lambda1p01_m4p5e04'):
    output_folder = 'bluemoon_H30_lambda1p01_m4p5e04'

elif (spacecraft.lander_type == 'bluemoon_H30_lambda1p5_m1p5e04'):
    output_folder = 'bluemoon_H30_lambda1p5_m1p5e04'
elif (spacecraft.lander_type == 'bluemoon_H30_lambda1p5_m3e04'):
    output_folder = 'bluemoon_H30_lambda1p5_m3e04'
elif (spacecraft.lander_type == 'bluemoon_H30_lambda1p5_m4p5e04'):
    output_folder = 'bluemoon_H30_lambda1p5_m4p5e04'

elif (spacecraft.lander_type == 'bluemoon_H30_lambda2_m1p5e04'):
    output_folder = 'bluemoon_H30_lambda2_m1p5e04'
elif (spacecraft.lander_type == 'bluemoon_H30_lambda2_m3e04'):
    output_folder = 'bluemoon_H30_lambda2_m3e04'
elif (spacecraft.lander_type == 'bluemoon_H30_lambda2_m4p5e04'):
    output_folder = 'bluemoon_H30_lambda2_m4p5e04'

elif (spacecraft.lander_type == 'bluemoon_H30_lambda2p5_m1p5e04'):
    output_folder = 'bluemoon_H30_lambda2p5_m1p5e04'
elif (spacecraft.lander_type == 'bluemoon_H30_lambda2p5_m3e04'):
    output_folder = 'bluemoon_H30_lambda2p5_m3e04'
elif (spacecraft.lander_type == 'bluemoon_H30_lambda2p5_m4p5e04'):
    output_folder = 'bluemoon_H30_lambda2p5_m4p5e04'

elif (spacecraft.lander_type == 'bluemoon_H30_lambda3_m1p5e04'):
    output_folder = 'bluemoon_H30_lambda3_m1p5e04'
elif (spacecraft.lander_type == 'bluemoon_H30_lambda3_m3e04'):
    output_folder = 'bluemoon_H30_lambda3_m3e04'
elif (spacecraft.lander_type == 'bluemoon_H30_lambda3_m4p5e04'):
    output_folder = 'bluemoon_H30_lambda3_m4p5e04'

elif (spacecraft.lander_type == 'bluemoon_H30_lambda4_m1p5e04'):
    output_folder = 'bluemoon_H30_lambda4_m1p5e04'
elif (spacecraft.lander_type == 'bluemoon_H30_lambda4_m3e04'):
    output_folder = 'bluemoon_H30_lambda4_m3e04'
elif (spacecraft.lander_type == 'bluemoon_H30_lambda4_m4p5e04'):
    output_folder = 'bluemoon_H30_lambda4_m4p5e04'

elif (spacecraft.lander_type == 'bluemoon_H30_lambda5_m1p5e04'):
    output_folder = 'bluemoon_H30_lambda5_m1p5e04'
elif (spacecraft.lander_type == 'bluemoon_H30_lambda5_m3e04'):
    output_folder = 'bluemoon_H30_lambda5_m3e04'
elif (spacecraft.lander_type == 'bluemoon_H30_lambda5_m4p5e04'):
    output_folder = 'bluemoon_H30_lambda5_m4p5e04'

elif (spacecraft.lander_type == 'bluemoon_H30_lambda6_m1p5e04'):
    output_folder = 'bluemoon_H30_lambda6_m1p5e04'
elif (spacecraft.lander_type == 'bluemoon_H30_lambda6_m3e04'):
    output_folder = 'bluemoon_H30_lambda6_m3e04'
elif (spacecraft.lander_type == 'bluemoon_H30_lambda6_m4p5e04'):
    output_folder = 'bluemoon_H30_lambda6_m4p5e04'

elif (spacecraft.lander_type == 'bluemoon_H60_lambda1p01_m1p5e04'):
    output_folder = 'bluemoon_H60_lambda1p01_m1p5e04'
elif (spacecraft.lander_type == 'bluemoon_H60_lambda1p01_m3e04'):
    output_folder = 'bluemoon_H60_lambda1p01_m3e04'
elif (spacecraft.lander_type == 'bluemoon_H60_lambda1p01_m4p5e04'):
    output_folder = 'bluemoon_H60_lambda1p01_m4p5e04'

elif (spacecraft.lander_type == 'bluemoon_H60_lambda1p5_m1p5e04'):
    output_folder = 'bluemoon_H60_lambda1p5_m1p5e04'
elif (spacecraft.lander_type == 'bluemoon_H60_lambda1p5_m3e04'):
    output_folder = 'bluemoon_H60_lambda1p5_m3e04'
elif (spacecraft.lander_type == 'bluemoon_H60_lambda1p5_m4p5e04'):
    output_folder = 'bluemoon_H60_lambda1p5_m4p5e04'

elif (spacecraft.lander_type == 'bluemoon_H60_lambda2_m1p5e04'):
    output_folder = 'bluemoon_H60_lambda2_m1p5e04'
elif (spacecraft.lander_type == 'bluemoon_H60_lambda2_m3e04'):
    output_folder = 'bluemoon_H60_lambda2_m3e04'
elif (spacecraft.lander_type == 'bluemoon_H60_lambda2_m4p5e04'):
    output_folder = 'bluemoon_H60_lambda2_m4p5e04'

elif (spacecraft.lander_type == 'bluemoon_H60_lambda2p5_m1p5e04'):
    output_folder = 'bluemoon_H60_lambda2p5_m1p5e04'
elif (spacecraft.lander_type == 'bluemoon_H60_lambda2p5_m3e04'):
    output_folder = 'bluemoon_H60_lambda2p5_m3e04'
elif (spacecraft.lander_type == 'bluemoon_H60_lambda2p5_m4p5e04'):
    output_folder = 'bluemoon_H60_lambda2p5_m4p5e04'

elif (spacecraft.lander_type == 'bluemoon_H60_lambda3_m1p5e04'):
    output_folder = 'bluemoon_H60_lambda3_m1p5e04'
elif (spacecraft.lander_type == 'bluemoon_H60_lambda3_m3e04'):
    output_folder = 'bluemoon_H60_lambda3_m3e04'
elif (spacecraft.lander_type == 'bluemoon_H60_lambda3_m4p5e04'):
    output_folder = 'bluemoon_H60_lambda3_m4p5e04'

elif (spacecraft.lander_type == 'bluemoon_H60_lambda4_m1p5e04'):
    output_folder = 'bluemoon_H60_lambda4_m1p5e04'
elif (spacecraft.lander_type == 'bluemoon_H60_lambda4_m3e04'):
    output_folder = 'bluemoon_H60_lambda4_m3e04'
elif (spacecraft.lander_type == 'bluemoon_H60_lambda4_m4p5e04'):
    output_folder = 'bluemoon_H60_lambda4_m4p5e04'

elif (spacecraft.lander_type == 'bluemoon_H60_lambda5_m1p5e04'):
    output_folder = 'bluemoon_H60_lambda5_m1p5e04'
elif (spacecraft.lander_type == 'bluemoon_H60_lambda5_m3e04'):
    output_folder = 'bluemoon_H60_lambda5_m3e04'
elif (spacecraft.lander_type == 'bluemoon_H60_lambda5_m4p5e04'):
    output_folder = 'bluemoon_H60_lambda5_m4p5e04'

elif (spacecraft.lander_type == 'bluemoon_H60_lambda6_m1p5e04'):
    output_folder = 'bluemoon_H60_lambda6_m1p5e04'
elif (spacecraft.lander_type == 'bluemoon_H60_lambda6_m3e04'):
    output_folder = 'bluemoon_H60_lambda6_m3e04'
elif (spacecraft.lander_type == 'bluemoon_H60_lambda6_m4p5e04'):
    output_folder = 'bluemoon_H60_lambda6_m4p5e04'

elif (spacecraft.lander_type == 'starship_H30_lambda1p5_m5e04'):
    output_folder = 'starship_H30_lambda1p5_m5e04'
elif (spacecraft.lander_type == 'starship_H30_lambda2_m5e04'):
    output_folder = 'starship_H30_lambda2_m5e04'
elif (spacecraft.lander_type == 'starship_H30_lambda4_m5e04'):
    output_folder = 'starship_H30_lambda4_m5e04'
elif (spacecraft.lander_type == 'starship_H30_lambda6_m5e04'):
    output_folder = 'starship_H30_lambda6_m5e04'
elif (spacecraft.lander_type == 'starship_H30_lambda8_m5e04'):
    output_folder = 'starship_H30_lambda8_m5e04'
elif (spacecraft.lander_type == 'starship_H30_lambda10_m5e04'):
    output_folder = 'starship_H30_lambda10_m5e04'

# -----------------------------------------------------------------
# ADDED Nov 1st 2025
# -----------------------------------------------------------------

elif (spacecraft.lander_type == 'bluemoon_H30_lambda1p5_m2p25e04'):
    output_folder = 'bluemoon_H30_lambda1p5_m2p25e04'
elif (spacecraft.lander_type == 'bluemoon_H30_lambda1p5_m3p75e04'):
    output_folder = 'bluemoon_H30_lambda1p5_m3p75e04'
elif (spacecraft.lander_type == 'bluemoon_H60_lambda1p5_m2p25e04'):
    output_folder = 'bluemoon_H60_lambda1p5_m2p25e04'
elif (spacecraft.lander_type == 'bluemoon_H60_lambda1p5_m3p75e04'):
    output_folder = 'bluemoon_H60_lambda1p5_m3p75e04'

elif (spacecraft.lander_type == 'bluemoon_H45_lambda1p5_m1p5e04'):
    output_folder = 'bluemoon_H45_lambda1p5_m1p5e04'
elif (spacecraft.lander_type == 'bluemoon_H45_lambda1p5_m2p25e04'):
    output_folder = 'bluemoon_H45_lambda1p5_m2p25e04'
elif (spacecraft.lander_type == 'bluemoon_H45_lambda1p5_m3e04'):
    output_folder = 'bluemoon_H45_lambda1p5_m3e04'
elif (spacecraft.lander_type == 'bluemoon_H45_lambda1p5_m3p75e04'):
    output_folder = 'bluemoon_H45_lambda1p5_m3p75e04'
elif (spacecraft.lander_type == 'bluemoon_H45_lambda1p5_m4p5e04'):
    output_folder = 'bluemoon_H45_lambda1p5_m4p5e04'
elif (spacecraft.lander_type == 'bluemoon_H90_lambda1p5_m1p5e04'):
    output_folder = 'bluemoon_H90_lambda1p5_m1p5e04'
elif (spacecraft.lander_type == 'bluemoon_H90_lambda1p5_m2p25e04'):
    output_folder = 'bluemoon_H90_lambda1p5_m2p25e04'
elif (spacecraft.lander_type == 'bluemoon_H90_lambda1p5_m3e04'):
    output_folder = 'bluemoon_H90_lambda1p5_m3e04'
elif (spacecraft.lander_type == 'bluemoon_H90_lambda1p5_m3p75e04'):
    output_folder = 'bluemoon_H90_lambda1p5_m3p75e04'
elif (spacecraft.lander_type == 'bluemoon_H90_lambda1p5_m4p5e04'):
    output_folder = 'bluemoon_H90_lambda1p5_m4p5e04'

elif (spacecraft.lander_type == 'bluemoon_H45_lambda2_m1p5e04'):
    output_folder = 'bluemoon_H45_lambda2_m1p5e04'
elif (spacecraft.lander_type == 'bluemoon_H45_lambda2p5_m1p5e04'):
    output_folder = 'bluemoon_H45_lambda2p5_m1p5e04'
elif (spacecraft.lander_type == 'bluemoon_H45_lambda3_m1p5e04'):
    output_folder = 'bluemoon_H45_lambda3_m1p5e04'
elif (spacecraft.lander_type == 'bluemoon_H45_lambda4_m1p5e04'):
    output_folder = 'bluemoon_H45_lambda4_m1p5e04'
elif (spacecraft.lander_type == 'bluemoon_H45_lambda5_m1p5e04'):
    output_folder = 'bluemoon_H45_lambda5_m1p5e04'
elif (spacecraft.lander_type == 'bluemoon_H90_lambda2_m1p5e04'):
    output_folder = 'bluemoon_H90_lambda2_m1p5e04'
elif (spacecraft.lander_type == 'bluemoon_H90_lambda2p5_m1p5e04'):
    output_folder = 'bluemoon_H90_lambda2p5_m1p5e04'
elif (spacecraft.lander_type == 'bluemoon_H90_lambda3_m1p5e04'):
    output_folder = 'bluemoon_H90_lambda3_m1p5e04'
elif (spacecraft.lander_type == 'bluemoon_H90_lambda4_m1p5e04'):
    output_folder = 'bluemoon_H90_lambda4_m1p5e04'
elif (spacecraft.lander_type == 'bluemoon_H90_lambda5_m1p5e04'):
    output_folder = 'bluemoon_H90_lambda5_m1p5e04'

# -----------------------------------------------------------------
# -----------------------------------------------------------------


elif (spacecraft.lander_type == 'starship_H60_lambda1p5_m5e04'):
    output_folder = 'starship_H60_lambda1p5_m5e04'
elif (spacecraft.lander_type == 'starship_H60_lambda2_m5e04'):
    output_folder = 'starship_H60_lambda2_m5e04'
elif (spacecraft.lander_type == 'starship_H60_lambda4_m5e04'):
    output_folder = 'starship_H60_lambda4_m5e04'
elif (spacecraft.lander_type == 'starship_H60_lambda6_m5e04'):
    output_folder = 'starship_H60_lambda6_m5e04'
elif (spacecraft.lander_type == 'starship_H60_lambda8_m5e04'):
    output_folder = 'starship_H60_lambda8_m5e04'
elif (spacecraft.lander_type == 'starship_H60_lambda10_m5e04'):
    output_folder = 'starship_H60_lambda10_m5e04'


elif (spacecraft.lander_type == 'starship_H30_lambda1p5_m1e05'):
    output_folder = 'starship_H30_lambda1p5_m1e05'
elif (spacecraft.lander_type == 'starship_H30_lambda2_m1e05'):
    output_folder = 'starship_H30_lambda2_m1e05'
elif (spacecraft.lander_type == 'starship_H30_lambda4_m1e05'):
    output_folder = 'starship_H30_lambda4_m1e05'
elif (spacecraft.lander_type == 'starship_H30_lambda6_m1e05'):
    output_folder = 'starship_H30_lambda6_m1e05'
elif (spacecraft.lander_type == 'starship_H30_lambda8_m1e05'):
    output_folder = 'starship_H30_lambda8_m1e05'
elif (spacecraft.lander_type == 'starship_H30_lambda10_m1e05'):
    output_folder = 'starship_H30_lambda10_m1e05'



elif (spacecraft.lander_type == 'starship_H60_lambda1p5_m1e05'):
    output_folder = 'starship_H60_lambda1p5_m1e05'
elif (spacecraft.lander_type == 'starship_H60_lambda2_m1e05'):
    output_folder = 'starship_H60_lambda2_m1e05'
elif (spacecraft.lander_type == 'starship_H60_lambda4_m1e05'):
    output_folder = 'starship_H60_lambda4_m1e05'
elif (spacecraft.lander_type == 'starship_H60_lambda6_m1e05'):
    output_folder = 'starship_H60_lambda6_m1e05'
elif (spacecraft.lander_type == 'starship_H60_lambda8_m1e05'):
    output_folder = 'starship_H60_lambda8_m1e05'
elif (spacecraft.lander_type == 'starship_H60_lambda10_m1e05'):
    output_folder = 'starship_H60_lambda10_m1e05'

# ----------------------------------------------------------------
# Starship Upper Thruster Routines
# ----------------------------------------------------------------

elif (spacecraft.lander_type == 'starshipUpper_H30_lambda1p5_m5e04'):
    output_folder = 'starshipUpper_H30_lambda1p5_m5e04'
elif (spacecraft.lander_type == 'starshipUpper_H30_lambda2_m5e04'):
    output_folder = 'starshipUpper_H30_lambda2_m5e04'
elif (spacecraft.lander_type == 'starshipUpper_H30_lambda4_m5e04'):
    output_folder = 'starshipUpper_H30_lambda4_m5e04'
elif (spacecraft.lander_type == 'starshipUpper_H30_lambda6_m5e04'):
    output_folder = 'starshipUpper_H30_lambda6_m5e04'
elif (spacecraft.lander_type == 'starshipUpper_H30_lambda8_m5e04'):
    output_folder = 'starshipUpper_H30_lambda8_m5e04'
elif (spacecraft.lander_type == 'starshipUpper_H30_lambda10_m5e04'):
    output_folder = 'starshipUpper_H30_lambda10_m5e04'

elif (spacecraft.lander_type == 'starshipUpper_H60_lambda1p5_m5e04'):
    output_folder = 'starshipUpper_H60_lambda1p5_m5e04'
elif (spacecraft.lander_type == 'starshipUpper_H60_lambda2_m5e04'):
    output_folder = 'starshipUpper_H60_lambda2_m5e04'
elif (spacecraft.lander_type == 'starshipUpper_H60_lambda4_m5e04'):
    output_folder = 'starshipUpper_H60_lambda4_m5e04'
elif (spacecraft.lander_type == 'starshipUpper_H60_lambda6_m5e04'):
    output_folder = 'starshipUpper_H60_lambda6_m5e04'
elif (spacecraft.lander_type == 'starshipUpper_H60_lambda8_m5e04'):
    output_folder = 'starshipUpper_H60_lambda8_m5e04'
elif (spacecraft.lander_type == 'starshipUpper_H60_lambda10_m5e04'):
    output_folder = 'starshipUpper_H60_lambda10_m5e04'

elif (spacecraft.lander_type == 'starshipUpper_H30_lambda1p5_m1e05'):
    output_folder = 'starshipUpper_H30_lambda1p5_m1e05'
elif (spacecraft.lander_type == 'starshipUpper_H30_lambda2_m1e05'):
    output_folder = 'starshipUpper_H30_lambda2_m1e05'
elif (spacecraft.lander_type == 'starshipUpper_H30_lambda4_m1e05'):
    output_folder = 'starshipUpper_H30_lambda4_m1e05'
elif (spacecraft.lander_type == 'starshipUpper_H30_lambda6_m1e05'):
    output_folder = 'starshipUpper_H30_lambda6_m1e05'
elif (spacecraft.lander_type == 'starshipUpper_H30_lambda8_m1e05'):
    output_folder = 'starshipUpper_H30_lambda8_m1e05'
elif (spacecraft.lander_type == 'starshipUpper_H30_lambda10_m1e05'):
    output_folder = 'starshipUpper_H30_lambda10_m1e05'

elif (spacecraft.lander_type == 'starshipUpper_H60_lambda1p5_m1e05'):
    output_folder = 'starshipUpper_H60_lambda1p5_m1e05'
elif (spacecraft.lander_type == 'starshipUpper_H60_lambda2_m1e05'):
    output_folder = 'starshipUpper_H60_lambda2_m1e05'
elif (spacecraft.lander_type == 'starshipUpper_H60_lambda4_m1e05'):
    output_folder = 'starshipUpper_H60_lambda4_m1e05'
elif (spacecraft.lander_type == 'starshipUpper_H60_lambda6_m1e05'):
    output_folder = 'starshipUpper_H60_lambda6_m1e05'
elif (spacecraft.lander_type == 'starshipUpper_H60_lambda8_m1e05'):
    output_folder = 'starshipUpper_H60_lambda8_m1e05'
elif (spacecraft.lander_type == 'starshipUpper_H60_lambda10_m1e05'):
    output_folder = 'starshipUpper_H60_lambda10_m1e05'

#if(spacecraft.lander_type == 'bluemoon_case_1'):
#    output_folder = 'bluemoon_case_1'
#if(spacecraft.lander_type == 'starship_raptor_nominal_50'):
#    output_folder = 'starship_raptor_nominal_50_scaled_dec2024'
#elif(spacecraft.lander_type == 'starship_raptor_nominal_100'):
#    output_folder = 'starship_raptor_nominal_100_scaled_dec2024' 
#elif(spacecraft.lander_type == 'starship_thruster_nominal_50'):
#    output_folder = 'starship_thruster_nominal_50_scaled'
#elif(spacecraft.lander_type == 'starship_thruster_nominal_100'):
#    output_folder = 'starship_thruster_nominal_100_scaled' 
#
#elif(spacecraft.lander_type == "apollo_vel_slow_alt_high"):
#    output_folder = "apollo_vel_slow_alt_high_4"
#elif(spacecraft.lander_type == "apollo_vel_slow_alt_mid"):
#    output_folder = "apollo_vel_slow_alt_mid_4"
#elif(spacecraft.lander_type == "apollo_vel_slow_alt_low"):
#    output_folder = "apollo_vel_slow_alt_low_4"
#elif(spacecraft.lander_type == "apollo_vel_mid_alt_mid"):
#    output_folder = "apollo_vel_mid_alt_mid_4"
#elif(spacecraft.lander_type == "apollo_vel_fast_alt_mid"):
#    output_folder = "apollo_vel_fast_alt_mid_4"
#
#elif(spacecraft.lander_type == 'starship_thruster_vel_high_100'):
#    output_folder = 'starship_thruster_vel_high_100' 
#elif(spacecraft.lander_type == 'starship_thruster_vel_high_50'):
#    output_folder = 'starship_thruster_vel_high_50' 
#elif(spacecraft.lander_type == 'starship_thruster_vel_mid_100'):
#    output_folder = 'starship_thruster_vel_mid_100' 
#elif(spacecraft.lander_type == 'starship_thruster_vel_low_100'):
#    output_folder = 'starship_thruster_vel_low_100' 
#
#
#elif(spacecraft.lander_type == 'starship_raptor_vel_high_100'):
#    output_folder = 'starship_raptor_vel_high_100' 
#elif(spacecraft.lander_type == 'starship_raptor_vel_high_50'):
#    output_folder = 'starship_raptor_vel_high_50' 
#elif(spacecraft.lander_type == 'starship_raptor_vel_mid_100'):
#    output_folder = 'starship_raptor_vel_mid_100' 
#elif(spacecraft.lander_type == 'starship_raptor_vel_low_100'):
#    output_folder = 'starship_raptor_vel_low_100' 

else:
    output_folder = 'output/'
    
output_folder = output_folder + "_scaling_"+str(scale.scaling_index)+"/"


def print_timestep(ts, h_nozzle):

    Nt = timestep.n_sub_timesteps

    if (ts < 10) or (ts % 10 == 0) or (ts >= Nt - 10):
        
        # add a blank line for spacing
        print()
        print("--------------------------------------------------------------------------------------------")
    
        # print the timestep update on one line
        print(
            f"Timestep {ts}/{Nt}  "
            f"t = {ts * timestep.delta_t:.3f} s  "
            f"Height = {h_nozzle:.3f} m  "
            f"Descent = {nozzle.v_descent:.5f} m/s  "
            f"S = {scale.scaling_factor[scale.scaling_index]}",
            flush=True
        )
        print("--------------------------------------------------------------------------------------------")
    

def plot_erosion_profile(ts, r_from_centerline, d_excavated):

    if not out.plot_profile:
        return None

    Nt = timestep.n_sub_timesteps

    # plot first 10, every 10, and last 10
    should_plot = (ts < 10) or (ts % 10 == 0) or (ts >= Nt - 10)

    if not should_plot:
        return None

    # clean line before plot print (so it doesn't overwrite timestep print)
    print()

    print(f"Plotting Figure at Timestep: {ts}")

    fig = plt.figure(figsize=(10, 6))
    plt.title(f"Time = {ts * timestep.delta_t:.2f} s", fontsize=18)
    plt.scatter(r_from_centerline, d_excavated)

    plt.xlabel("Distance from Plume Centerline (m)", fontsize=18)
    plt.ylabel("Excavation Depth (m)", fontsize=18)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)

    plt.xlim(0, bounds.max_centerline)
    plt.ylim(0.5, 0)

    filename = f"{output_folder}{spacecraft.lander_type}_growth_{ts}.png"
    plt.savefig(filename, dpi=300)
    plt.close()

    return None


def write_bounds(bounds_ej_ring):

    if(out.write_bounds):
        np.savetxt(output_folder+spacecraft.lander_type+"_ejecta_ring_bounds.csv", bounds_ej_ring, delimiter=',', fmt=' '.join(['%i'] + ['%.4e']*3), header = '# Ejecta Ring Radial Distances Away from Plume Centerline (Index, Min (m), Midpoint (m), Max (m))')
        #print("Timestep,", "Depth Excavated,", "Threshold Energy,", "E_down,", "Alpha,","Mdot_flux", "M_area_eroded_inst", "Mdot_cumulative" )

    return None


def write_ejecta_props(ts, h_nozzle, ej_timestep_props, u_ej_timestep, offset_ej_dist_timestep, offset_ej_time_timestep):
    
    # saving values
    if(out.write_props):
        np.savetxt(output_folder+spacecraft.lander_type+"_ejecta_properties_"+str(ts)+".csv", ej_timestep_props, delimiter=',', fmt=' '.join(['%i'] + ['%.8e']*4), header = '# Nozzle Height (m) ' + str(h_nozzle)+ ', Descent Velocity (m/s) '+str(nozzle.v_descent)+', Offset angle (3 deg), ' + 'Bin Index, Bin Midpoint Location (m), Excavated Density (kg/m^3), Instantaneous Mass Excavated (kg), Total Height Excavated (m)')
        np.savetxt(output_folder+spacecraft.lander_type+"_ejecta_velocities_"+str(ts)+".csv", u_ej_timestep, delimiter=',', fmt=' '.join(['%i'] + ['%.4e']*np.size(soil.d_particle)), header = '# Nozzle Height (m) ' + str(h_nozzle)+ ', Descent Velocity (m/s) '+str(nozzle.v_descent)+', Offset angle (3 deg), '+' Ejecta Speeds by Bin Index and Particle Size (O(1), O(10), O(100), O(1000) microns) in sets of 1, 2, 3, 4, 5, 6, 7, 8, 9')
        np.savetxt(output_folder+spacecraft.lander_type+"_ejecta_offset_distances_"+str(ts)+".csv", offset_ej_dist_timestep, delimiter=',', fmt=' '.join(['%i'] + ['%.4e']*np.size(soil.d_particle)), header = '# header = Nozzle Height (m) ' + str(h_nozzle)+', Descent Velocity (m/s) '+str(nozzle.v_descent)+', Offset angle (3 deg), ' + 'Ejecta Distance Offset by Bin Index and Particle Size (O(1), O(10), O(100), O(1000) microns) in sets of 1, 2, 3, 4, 5, 6, 7, 8, 9')
        np.savetxt(output_folder+spacecraft.lander_type+"_ejecta_offset_time_"+str(ts)+".csv", offset_ej_time_timestep, delimiter=',', fmt=' '.join(['%i'] + ['%.4e']*np.size(soil.d_particle)), header = '# header = Nozzle Height (m) ' + str(h_nozzle)+', Descent Velocity (m/s) '+str(nozzle.v_descent)+', Offset angle (3 deg), ' + 'Ejecta Time Offset by Bin Index and Particle Size (O(1), O(10), O(100), O(1000) microns) in sets of 1, 2, 3, 4, 5, 6, 7, 8, 9')

    return None

def output_gas_props():
    if(out.output_gas):
        for i in range(0,bounds.n_points_centerline-1):
            print(i, soil.r_midpoint[i], imp.v_gas_arr[i], imp.p_gas_arr[i], imp.rho_gas_arr[i], imp.T_gas_arr[i])

    return None

def print_simulation_setup():
    """
    Prints all key input and derived parameters for the nozzle and timestepping setup.
    Thrust is computed assuming ambient pressure (P_amb) = 0 (vacuum).
    """

    pamb = 0.0  # Lunar vacuum ambient pressure [Pa]

    print("\n================= Simulation Setup =================")
    print(f"Lander type:                {nozzle.spacecraft.lander_type}")
    print(f"Engines (n):                {nozzle.n_engines}")
    print(f"Total Lander Mass (kg):     {nozzle.m_lander}")

    print("\nGeometry")
    print(f"  A_nozzle      [m^2]:      {nozzle.A_nozzle:.6g}")
    print(f"  A_throat      [m^2]:      {nozzle.A_throat:.6g}")
    print(f"  D_nozzle        [m]:      {nozzle.D_nozzle:.6g}")
    print(f"  r_nozzle        [m]:      {nozzle.r_nozzle:.6g}")

    print("\nKinematics / Descent")
    print(f"  h_nozzle_init   [m]:      {nozzle.h_nozzle_init:.6g}")
    print(f"  min_altitude    [m]:      {nozzle.min_altitude_lander:.6g}")
    print(f"  v_init        [m/s]:      {nozzle.v_init:.6g}  (downward)")
    print(f"  a_thrust    [m/s^2]:      {nozzle.a_thrust:.6g} (upward)")
    print(f"  total_sim_time   [s]:     {nozzle.lander_total_sim_time:.6g}")

    print("\nChamber Conditions")
    print(f"  P0             [Pa]:      {nozzle.P_0:.6g}")
    print(f"  T0              [K]:      {nozzle.T_0:.6g}")
    print(f"  rho0      [kg/m^3]:      {nozzle.rho_0:.6g}")

    print("\nExit Conditions")
    print(f"  Pe             [Pa]:      {nozzle.P_e:.6g}")
    print(f"  Te              [K]:      {nozzle.T_e:.6g}")
    print(f"  rho_e     [kg/m^3]:      {nozzle.rho_e:.6g}")
    print(f"  m_dot_e     [kg/s]:      {nozzle.m_dot_e:.6g}")
    print(f"  gamma           [-]:      {nozzle.gamma:.6g}")
    print(f"  Ma              [-]:      {nozzle.Ma:.6g}")

    print("\nDerived Exit Vars")
    print(f"  R_gas     [J/(kg·K)]:     {nozzle.R_gas:.6g}")
    print(f"  v_e            [m/s]:     {nozzle.v_e:.6g}")
    print(f"  k_bar           [-]:      {nozzle.k_bar:.6g}")

    # Compute thrust
    thrust_per_engine = (
        (nozzle.m_dot_e / nozzle.n_engines) * nozzle.v_e
        + (nozzle.P_e - pamb) * (nozzle.A_nozzle / nozzle.n_engines)
    )
    total_thrust = (
        nozzle.m_dot_e * nozzle.v_e
        + (nozzle.P_e - pamb) * nozzle.A_nozzle
    )

    lander_grav_force = nozzle.m_lander * constants.g 


    print("\nThrust Forces (Vacuum, P_amb = 0)")
    print(f"  Thrust per engine [N]:    {thrust_per_engine:.6g}")
    print(f"  Total thrust      [N]:    {total_thrust:.6g}")
    print(f"  Lander Grav Force [N]:    {lander_grav_force:.6g}")

    print("\nTimestepping")
    print(f"  Δt main        [s]:       {timestep.delta_t}")
    print(f"  dt substep     [s]:       {timestep.dt}")
    print(f"  n_sub_timesteps [-]:      {timestep.n_sub_timesteps}")
    print("====================================================\n")

#def plot_deposited_props():
#
#
#    if(sol_ej.plot_only):
#        sol_ej.range_bounds_mid = np.genfromtxt("output/range_bounds_"+str(sol_ej.index_grain_sizes[0])+"_"+str(sol_ej.index_grain_sizes[1])+".csv", delimiter=",")
#        sol_ej.mass_flux_upon_impact = np.genfromtxt("output/mass_flux_upon_impact_"+str(sol_ej.index_grain_sizes[0])+"_"+str(sol_ej.index_grain_sizes[1])+".csv", delimiter=',')
#        sol_ej.momentum_flux_upon_impact = np.genfromtxt("output/momentum_flux_upon_impact_"+str(sol_ej.index_grain_sizes[0])+"_"+str(sol_ej.index_grain_sizes[1])+".csv", delimiter=',')
#        sol_ej.energy_flux_upon_impact = np.genfromtxt("output/energy_flux_upon_impact_"+str(sol_ej.index_grain_sizes[0])+"_"+str(sol_ej.index_grain_sizes[1])+".csv", delimiter=',')
#        sol_ej.count_flux_upon_impact = np.genfromtxt("output/count_flux_upon_impact_"+str(sol_ej.index_grain_sizes[0])+"_"+str(sol_ej.index_grain_sizes[1])+".csv", delimiter=',')
#
#    total_mass_flux_upon_impact = np.zeros(sol_ej.n_sorted_bins)
#    total_momentum_flux_upon_impact = np.zeros(sol_ej.n_sorted_bins)
#    total_energy_flux_upon_impact = np.zeros(sol_ej.n_sorted_bins)
#    total_count_flux_upon_impact = np.zeros(sol_ej.n_sorted_bins)
#
#    for i in range(0,sol_ej.n_sorted_bins):
#        total_mass_flux_upon_impact[i] = np.sum(sol_ej.mass_flux_upon_impact[:,i])
#        total_momentum_flux_upon_impact[i] = np.sum(sol_ej.momentum_flux_upon_impact[:,i])
#        total_energy_flux_upon_impact[i] = np.sum(sol_ej.energy_flux_upon_impact[:,i])
#        total_count_flux_upon_impact[i] = np.sum(sol_ej.count_flux_upon_impact[:,i])
#
#    plt.figure(figsize=(5,5), layout = 'constrained')
#    plt.plot(sol_ej.range_bounds_mid, total_energy_flux_upon_impact, color = 'black', linewidth = 2, label = 'Total')
#
#    for i in range(sol_ej.index_grain_sizes[0],sol_ej.index_grain_sizes[1]):
#        #plt.ylim(0,10)
#        plt.plot(sol_ej.range_bounds_mid, sol_ej.energy_flux_upon_impact[i,:], linewidth = 2, linestyle = 'dotted', label = "Particle Size (m): " + str('{:.1E}'.format(soil.d_particle[i])))
#        #plt.xscale('log')
#        plt.yscale('log')
#
#    plt.xlabel("Distance from Centerline (m)", fontsize = 16)
#    plt.ylabel("Impact Energy $(J/m^2)$", fontsize = 16)
#    plt.xticks(fontsize = 16)
#    plt.yticks(fontsize = 16)
#    plt.locator_params(axis='x', nbins=5)
#    plt.xlim(0,sol_ej.range_of_analysis)
#    #plt.ylim(1E-1,1E3)
#    plt.grid()
#    plt.legend(fontsize = 8)
#    plt.savefig("output/impact_energy_profile_"+str(sol_ej.index_grain_sizes[0])+"_"+str(sol_ej.index_grain_sizes[1])+".png", bbox_inches='tight', dpi=100)
#    plt.close()
#
#
#    plt.figure(figsize=(5,5), layout = 'constrained')
#    plt.plot(sol_ej.range_bounds_mid, total_momentum_flux_upon_impact, color = 'black', linewidth = 2, label = 'Total')
#
#    for i in range(sol_ej.index_grain_sizes[0], sol_ej.index_grain_sizes[1]):
#        #plt.ylim(0,10)
#        plt.plot(sol_ej.range_bounds_mid, sol_ej.momentum_flux_upon_impact[i,:], linewidth = 2, linestyle = 'dotted', label = "Particle Size (m): " + str('{:.1E}'.format(soil.d_particle[i])))
#        #plt.xscale('log')
#        plt.yscale('log')
#
#    plt.xlabel("Distance from Centerline (m)", fontsize = 16)
#    plt.ylabel("Impact Momentum $(N \cdot s/m^2)$", fontsize = 16)
#    plt.xticks(fontsize = 16)
#    plt.yticks(fontsize = 16)
#    plt.locator_params(axis='x', nbins=5)
#    plt.xlim(0,sol_ej.range_of_analysis)
#    plt.grid()
#    plt.legend(fontsize = 8)
#    plt.savefig("output/impact_momentum_profile_"+str(sol_ej.index_grain_sizes[0])+"_"+str(sol_ej.index_grain_sizes[1])+".png", bbox_inches='tight', dpi=100)
#    plt.close()
#
#    plt.figure(figsize=(5,5), layout = 'constrained')
#    plt.plot(sol_ej.range_bounds_mid, total_mass_flux_upon_impact, color = 'black', linewidth = 2, label = 'Total')
#
#    for i in range(sol_ej.index_grain_sizes[0],sol_ej.index_grain_sizes[1]):
#        #plt.ylim(0,10)
#        plt.plot(sol_ej.range_bounds_mid, sol_ej.mass_flux_upon_impact[i,:], linewidth = 2, linestyle = 'dotted', label = "Particle Size (m): " + str('{:.1E}'.format(soil.d_particle[i])))
#        #plt.xscale('log')
#        plt.yscale('log')
#
#    plt.xlabel("Distance from Centerline (m)", fontsize = 16)
#    plt.ylabel("Impact Mass $(kg/m^2)$", fontsize = 16)
#    plt.xticks(fontsize = 16)
#    plt.yticks(fontsize = 16)
#    plt.locator_params(axis='x', nbins=5)
#    plt.xlim(0,sol_ej.range_of_analysis)
#    plt.grid()
#    plt.legend(fontsize = 8)
#    plt.savefig("output/impact_mass_profile_"+str(sol_ej.index_grain_sizes[0])+"_"+str(sol_ej.index_grain_sizes[1])+".png", bbox_inches='tight',dpi=100)
#    plt.close()
#
#
#    plt.figure(figsize=(5,5), layout = 'constrained')
#    plt.plot(sol_ej.range_bounds_mid, total_count_flux_upon_impact, color = 'black', linewidth = 2, label = 'Total')
#
#    for i in range(sol_ej.index_grain_sizes[0],sol_ej.index_grain_sizes[1]):
#        #plt.ylim(0,10)
#        plt.plot(sol_ej.range_bounds_mid, sol_ej.count_flux_upon_impact[i,:], linewidth = 2, linestyle = 'dotted', label = "Particle Size (m): " + str('{:.1E}'.format(soil.d_particle[i])))
#        #plt.xscale('log')
#        plt.yscale('log')
#
#    plt.xlabel("Distance from Centerline (m)", fontsize = 16)
#    plt.ylabel("Impact Count $(\#/m^2)$", fontsize = 16)
#    plt.xticks(fontsize = 16)
#    plt.yticks(fontsize = 16)
#    plt.locator_params(axis='x', nbins=5)
#    plt.xlim(0,sol_ej.range_of_analysis)
#    plt.grid()
#    plt.legend(fontsize = 8)
#    plt.savefig("output/impact_count_profile_"+str(sol_ej.index_grain_sizes[0])+"_"+str(sol_ej.index_grain_sizes[1])+".png", bbox_inches='tight',dpi=100)
#    plt.close()
#
#
#    plt.figure(figsize=(5,5), layout = 'constrained')
#    plt.plot(sol_ej.range_bounds_mid, total_count_flux_upon_impact * 0.0025, color = 'black', linewidth = 2, label = 'Total')
#
#    for i in range(sol_ej.index_grain_sizes[0],sol_ej.index_grain_sizes[1]):
#        #plt.ylim(0,10)
#        plt.plot(sol_ej.range_bounds_mid, sol_ej.count_flux_upon_impact[i,:] * 0.0025, linewidth = 2, linestyle = 'dotted', label = "Particle Size (m): " + str('{:.1E}'.format(soil.d_particle[i])))
#        #plt.xscale('log')
#        plt.yscale('log')
#
#    plt.xlabel("Distance from Centerline (m)", fontsize = 16)
#    plt.ylabel("Impact Count $(\#/25cm^2)$", fontsize = 16)
#    plt.xticks(fontsize = 16)
#    plt.yticks(fontsize = 16)
#    plt.locator_params(axis='x', nbins=5)
#    plt.xlim(0,sol_ej.range_of_analysis)
#    plt.grid()
#    plt.legend(fontsize = 8)
#    plt.savefig("output/impact_count_sample_profile_"+str(sol_ej.index_grain_sizes[0])+"_"+str(sol_ej.index_grain_sizes[1])+".png", bbox_inches='tight',dpi=100)
#    plt.close()



    return None

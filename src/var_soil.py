# universal variables
#rho_p = 1800
        
# Calculate the surface threshold parameters
E_th_0 = 0.123 #J/m^2/s
epsilon = 0.0029 # average value for soil efficiency
#alpha = 0.289 # J/m^3 (constant)
alpha_0 = 0.289
D_84 = 450E-6 # m, Check back with apollo 16 data for more accurate answer
D_bracket = 1.5 * D_84

rho_0 = 1100
rho_inf = 1800 #38% porosity deep within crater
h = 8.0/100 # From Metzger (2024)

k_cohesion = 50 #kg/m^3 (can change)

rho_bulk = None
E_downward = None
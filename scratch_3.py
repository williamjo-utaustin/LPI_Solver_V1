import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

d_particle_array = np.asarray([multiplier * magnitude for magnitude in [1E-6, 1E-5, 1E-4, 1E-3] for multiplier in [1, 2, 3, 4, 5, 6, 7, 8, 9]])

def dudt(y, t, d_p, u_gas, p_gas, rho_gas, T_gas, rho_p):
    g = 1.62      # Gravity m/s^2
    A_const= 1.71575E-7
    beta = 0.78


    mfp = np.sqrt(np.pi/(2 * p_gas * rho_gas)) * A_const * T_gas**beta
    Kn_p = mfp/d_p
    Re = (rho_gas * np.abs(u_gas - y) * d_p)/(A_const* T_gas ** beta)

    print(Re)


    C_D = cd_sphere(Re, Kn_p)
    #F_D = (3.0*rho_gas*C_D)*(u_gas - y)**2/(4.0*rho_p*d_p)
    F_D_stokes = 6 * np.pi * A_const*T_gas**beta * (d_p/2) * (u_gas - y)**2

    #F_D = 0.5 * rho_gas * (u_gas - y)**2 * C_D * (np.pi * (d_p/2)**2)

    #print(F_D + ((1 - (rho_gas/rho_p))*g))

    return F_D_stokes + ((1 - (rho_gas/rho_p))*g)


#def cd_sphere(Re, Kn):
#    
#    if Re <= 0.0:
#        CD = 0.0
#    elif Re > 0.0 and Re <= 2:
#        CD = (24.0/Re)
#    elif Re >= 2 and Re < 500:
#        CD = 18.5 * Re **-0.6
#    else:
#        CD = 0.44
#    
#    S_correction = 1 + Kn * (2.514 + 0.8 * np.exp(-0.55/Kn))
#    CD = CD/S_correction
#    return CD

def cd_sphere(Re, Kn):
    
    if Re <= 0.0:
        CD = 0.0
    elif Re > 8.0e6:
        CD = 0.2
    elif Re > 0.0 and Re <= 1:
        CD = (24.0/Re)
    elif Re > 1 and Re <= 100.0:
        p = np.array([4.22, -14.05, 34.87, 0.658])
        CD = np.polyval(p, 1.0/Re) 
    elif Re > 100.0 and Re <= 1.0e4:
        p = np.array([-30.41, 43.72, -17.08, 2.41])
        CD = np.polyval(p, 1.0/np.log10(Re))
    elif Re > 1.0e4 and Re <= 3.35e5:
        p = np.array([-0.1584, 2.031, -8.472, 11.932])
        CD = np.polyval(p, np.log10(Re))
    elif Re > 3.35e5 and Re <= 5.0e5:
        x1 = np.log10(Re/4.5e5)
        CD = 91.08*x1**4 + 0.0764
    else:
        p = np.array([-0.06338, 1.1905, -7.332, 14.93])
        CD = np.polyval(p, np.log10(Re))

    # add in correction factor to account for rarefied effects
    S_correction = 1 + Kn * (2.514 + 0.8 * np.exp(-0.55/Kn))
    CD = CD/S_correction
    return CD



rho_p = 1100.35185
rho_gas = 6.929E-11
p_gas = 0.001376
T_gas = 4762
u_gas = 3665

for i in range(0,np.size(d_particle_array)):
    d_p = d_particle_array[i]
    args_gas = (d_p, u_gas, p_gas, rho_gas, T_gas, rho_p)
   
    sol = integrate.solve_ivp(dudt, [0, 10000], [0, 0], method='DOP853', args = args_gas, rtol = 1E-10)
    #sol = integrate.odeint(dudt, 0, np.linspace(0,2000,2000000), args = args_gas, rtol = 1E-11)


    #plt.plot(np.linspace(0,1000,2000000), sol[:,0])
    #plt.show()
    plt.plot(np.array(sol.t), np.array(sol.y[0,:]))
    plt.show()
    
    #plt.plot(np.array(sol.t), np.array(sol.y))
    #plt.show()

    
    



#x0 = 0.8
#args = (a,)
#y = integrate.odeint(func, x0, t, args)
#
#fig = plt.figure()
#ax = fig.add_subplot(111)
#h1, = ax.plot(t, y)
#h2, = ax.plot(t, [a(s) for s in t])
#ax.legend([h1, h2], ["y", "a"])
#ax.set_xlabel("t")
#ax.grid()
#plt.show()
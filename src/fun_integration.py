import numpy as np
import var_constants as cs

def rk_4(func, u_n, t, dt, d_p, p_gas, rho_gas, T_gas, u_gas, rho_p):
    dt2 = dt/2.0
    k1 = func(u_n, t, d_p, p_gas, rho_gas, T_gas, u_gas, rho_p)
    k2 = func(u_n + k1*dt2, t + dt2, d_p, p_gas, rho_gas, T_gas, u_gas, rho_p)
    k3 = func(u_n + k2*dt2, t + dt2, d_p, p_gas, rho_gas, T_gas, u_gas, rho_p)
    k4 = func(u_n + k3*dt, t + dt, d_p, p_gas, rho_gas, T_gas, u_gas, rho_p)
    u_n_p1 = u_n + dt/6.0*(k1 + 2.0*k2 + 2.0*k3 + k4)
    return u_n_p1

def dudt(u, t, d_p, p_gas, rho_gas, T_gas, u_gas, rho_p):
    A_const= 1.71575E-7
    beta = 0.78
    mean_free_path = np.sqrt(np.pi/(2 * p_gas * rho_gas)) * A_const * T_gas**beta
    Kn_p = mean_free_path/d_p
    Re = (rho_gas * np.abs(u_gas - u) * d_p)/(A_const* T_gas ** beta)
    C_D = cd_sphere(Re, Kn_p)
    F_D = (3.0*rho_gas*C_D)*(u_gas - u)**2/(4.0*rho_p*d_p)
    return F_D + ((1 - (rho_gas/rho_p))*cs.g)


def cd_sphere(Re, Kn):
    
    if Re <= 0.0:
        CD = 0.0
    elif Re > 8.0e6:
        CD = 0.2
    elif Re > 0.0 and Re <= 2:
        CD = (24.0/Re)
    elif Re > 2 and Re <= 100.0:
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
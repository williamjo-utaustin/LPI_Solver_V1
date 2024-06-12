import numpy as np

c = np.array([1.362, 0.415, 3.718, 5.084, 1.846, 0.770, 0.62, 0.4, 0.695])

param = np.asarray([multiplier * magnitude for magnitude in [1E0, 1E1, 1E2, 1E3] for multiplier in [1, 2, 3, 4, 5, 6, 7, 8, 9]])

def v_p(d_p, r_p, l_h, l_m):
    alpha = l_m**c[0] * l_h**(-c[1]) * d_p**(-c[2])
    #print(alpha)
    beta = c[3] * l_h * np.exp(-(r_p - (l_h - c[4]))**2/(2.0*c[5]**2))
    print(beta)
    gamma = np.exp(-(r_p - (l_h**c[6] - c[7]))**2/(2.0*c[8]**2))
    print(gamma)
    velocity = alpha * (beta + gamma)
    return velocity

u_p = v_p(param, 4000, 4000, 11000)
print(u_p)

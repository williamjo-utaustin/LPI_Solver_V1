import numpy as np

def A_over_Astar(M, gamma):
    return 1 / M * ((1 + (gamma - 1) / 2 * M ** 2) / ((gamma + 1) / 2)) ** ((gamma + 1) / (2 * (gamma - 1)))


def mDot_over_A(M, gamma, p0, R, T0):
    const = p0 / np.sqrt(R * T0)
    f = np.sqrt(gamma) * M / (1 + (gamma - 1) / 2 * M ** 2) ** ((gamma + 1) / (2 * (gamma - 1)))
    return const * f


def T0_over_T(M, gamma):
    return 1 + (gamma - 1) / 2 * M ** 2


def p0_over_p(M, gamma):
    return (T0_over_T(M, gamma)) ** (gamma / (gamma - 1))


def rho0_over_rho(M, gamma):
    return (T0_over_T(M, gamma)) ** (1 / (gamma - 1))

def pNS_over_p0(M, gamma):
    return (1 + gamma * M**2)/p0_over_p(M, gamma)


import numpy as np
from scipy.special import gamma
from scipy.stats import poisson, uniform, gamma as gamma_dist, moment


def theorCumulants(X0, alpha, beta, c, dt, b, flag):
    """
    Theoretical cumulants computation for CTS-OU and OU-CTS as described in: [3] Sabino, [4] Sabino

    INPUT
    X0:      initial condition
    alpha:   stability parameter
    beta:    beta
    c:       model parameter
    dt:      time interval
    b:       mean reverting parameter
    flag:    1 -> Finite Activity
             2 -> CTS-OU Finite Variation
             3 -> OU-CTS Finite Variation

    OUTPUT
    cumulants: computed cumulants as a numpy array
    """

    # Quantities of interest
    cumulants = np.zeros(4)
    k = np.arange(1, 5)

    # Cumulants computation
    if flag == 1:
        cumulants_L = c * beta**(alpha-k) * gamma(k-alpha)
        cumulants = X0 * np.exp(-b * dt) * (k == 1) + cumulants_L / (b * k) * (1 - np.exp(-k * b * dt))
    elif flag == 2:
        cumulants_X = c * beta ** (alpha - k) * gamma(k - alpha)
        cumulants = X0 * np.exp(-b * dt) * (k == 1) + cumulants_X * (1 - np.exp(-k * b * dt))
    elif flag == 3:
        cumulants_L = c * beta ** (alpha - k) * gamma(k - alpha)
        cumulants = X0 * np.exp(-b * dt) * (k == 1) + cumulants_L / (b * k) * (1 - np.exp(-k * b * dt))

    return cumulants

def bctsCumulants(X0, alpha, beta_p, beta_n, c_p, c_n, gamma_c, dt, b, flag):
    """
    Bilateral TS Cumulants Computation for CTS-OU and OU-CTS as described in: [1] Baviera

    INPUT
    X0:      initial condition
    alpha:   stability parameter
    beta_p:  positive beta
    beta_n:  negative beta
    c_p:     positive c
    c_n:     negative c
    gamma_c: drift
    dt:      time interval
    b:       mean reverting parameter
    flag:    1 -> OU-CTS Finite Activity
             2 -> CTS-OU Finite Variation
             3 -> OU-CTS Finite Variation

    OUTPUT
    cumulants: computed cumulants as a numpy array
    """

    # Quantities of interest
    cumulants = np.zeros(4)
    cumulants_L = np.zeros(4)
    cumulants_X = np.zeros(4)
    k = np.arange(1, 5)

    # Cumulants computation
    if flag == 1 or flag == 3:
        cumulants_L[0] = gamma_c
        cumulants_L[1:] = (c_p * beta_p**(alpha-k[1:]) * gamma(k[1:]-alpha) +
                           (-1)**k[1:] * c_n * beta_n**(alpha-k[1:]) * gamma(k[1:]-alpha))

        cumulants = (X0 * np.exp(-b * dt) * (k == 1) +
                     cumulants_L / (b * k) * (1 - np.exp(-k * b * dt)))
    elif flag == 2:
        cumulants_X[0] = gamma_c
        cumulants_X[1:] = (c_p * beta_p**(alpha-k[1:]) * gamma(k[1:]-alpha) +
                           (-1)**k[1:] * c_n * beta_n**(alpha-k[1:]) * gamma(k[1:]-alpha))

        cumulants = (X0 * np.exp(-b * dt) * (k == 1) +
                     cumulants_X * (1 - np.exp(-k * b * dt)))

    return cumulants

def compCumulants(vec):
    """
    Empirical cumulants computation

    INPUT
    vec: column vector of simulations

    OUTPUT
    cumulants: computed cumulants as a numpy array
    """
    # Quantities of interest
    cumulants = np.zeros(4)

    # Moments computation
    cumulants[0] = np.mean(vec)
    cumulants[1] = moment(vec, moment=2)
    cumulants[2] = moment(vec, moment=3)
    cumulants[3] = moment(vec, moment=4) - 3 * moment(vec, moment=2)**2

    return cumulants

def printCumulants(theorCumulants, simCumulants, T):
    print(f"\n Statistical cumulants with T {T:.2f} : ")
    print(f"Statistical cumulants order 1 = {simCumulants[0]:.5f}")
    print(f"Statistical cumulants order 2 = {simCumulants[1]:.5f}")
    print(f"Statistical cumulants order 3 = {simCumulants[2]:.5f}")
    print(f"Statistical cumulants order 4 = {simCumulants[3]:.5f}")
    print(f"\n Theoretical cumulants with delta t {T:.2f} : ")
    print(f"Theoretical cumulants order 1 = {theorCumulants[0]:.5f}")
    print(f"Theoretical cumulants order 2 = {theorCumulants[1]:.5f}")
    print(f"Theoretical cumulants order 3 = {theorCumulants[2]:.5f}")
    print(f"Theoretical cumulants order 4 = {theorCumulants[3]:.5f}")

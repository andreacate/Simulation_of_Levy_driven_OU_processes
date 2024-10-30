
import numpy as np
from scipy.integrate import quad
from scipy.special import gamma

def psi_X_TS_OU(u, alpha, beta_p, beta_n, c_p, c_n, gamma_c):
    """
    Characteristic function ([1] Baviera & Manzoni).

    Parameters:
    u       - Variable for which the characteristic function is computed (numpy array or scalar)
    alpha   - Shape parameter
    beta_p  - Positive beta parameter
    beta_n  - Negative beta parameter
    c_p     - Positive coefficient
    c_n     - Negative coefficient
    gamma_c - Constant shift parameter

    Returns:
    psi     - The computed characteristic function value
    """
    i = 1j

    term1 = i * u * gamma_c
    term2 = c_p * gamma(-alpha) * beta_p ** alpha * ((1 - i * u / beta_p) ** alpha - 1 + i * u * alpha / beta_p)
    term3 = c_n * gamma(-alpha) * beta_n ** alpha * ((1 + i * u / beta_n) ** alpha - 1 - i * u * alpha / beta_n)

    psi = term1 + term2 + term3
    return psi


def psi_Z_OU_TS(u, dt, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c):
    """
    Characteristic function ([1] Baviera & Manzoni)

    Parameters:
    u       - Variable for which the characteristic function is computed (numpy array or scalar)
    dt      - Time increment
    alpha   - Shape parameter
    b       - Mean reversion parameter
    beta_p  - Positive beta parameter
    beta_n  - Negative beta parameter
    c_p     - Positive coefficient
    c_n     - Negative coefficient
    gamma_c - Constant shift parameter

    Returns:
    psi     - The computed characteristic function value
    """

    u = np.atleast_1d(u)  # Ensure u is an array

    # Define integrands for real and imaginary parts
    def integrand_p_real(z, ui):
        return np.real((z - 1j * ui) ** alpha / z ** (alpha + 1))

    def integrand_p_imag(z, ui):
        return np.imag((z - 1j * ui) ** alpha / z ** (alpha + 1))

    def integrand_n_real(z, ui):
        return np.real((z + 1j * ui) ** alpha / z ** (alpha + 1))

    def integrand_n_imag(z, ui):
        return np.imag((z + 1j * ui) ** alpha / z ** (alpha + 1))

    # Calculate integrals for each value of u
    int_p_real = np.array(
        [quad(integrand_p_real, beta_p, beta_p * np.exp(b * dt), args = (ui,), limit = 100)[0] for ui in u])
    int_p_imag = np.array(
        [quad(integrand_p_imag, beta_p, beta_p * np.exp(b * dt), args = (ui,), limit = 100)[0] for ui in u])
    int_n_real = np.array(
        [quad(integrand_n_real, beta_n, beta_n * np.exp(b * dt), args = (ui,), limit = 100)[0] for ui in u])
    int_n_imag = np.array(
        [quad(integrand_n_imag, beta_n, beta_n * np.exp(b * dt), args = (ui,), limit = 100)[0] for ui in u])

    int_p = int_p_real + 1j * int_p_imag
    int_n = int_n_real + 1j * int_n_imag

    # Calculate characteristic exponent for each u
    term1 = 1j * u * (1 - np.exp(-b * dt)) / b * gamma_c
    term2 = c_p * beta_p ** alpha * gamma(-alpha) / b * (
                int_p - b * dt + alpha / beta_p * 1j * u * (1 - np.exp(-b * dt)))
    term3 = c_n * beta_n ** alpha * gamma(-alpha) / b * (
                int_n - b * dt - alpha / beta_n * 1j * u * (1 - np.exp(-b * dt)))

    psi = term1 + term2 + term3

    return psi


def psi_V(u, dt, alpha, b, beta_p, beta_n, c_p, c_n):
    """
       Characteristic function ([1] Baviera & Manzoni) present in the algorithm for the OU-TS finite activity

       Parameters:
       u       - Variable for which the characteristic function is computed (numpy array or scalar)
       dt      - Time increment
       alpha   - Shape parameter
       b       - Mean reversion parameter
       beta_p  - Positive beta parameter
       beta_n  - Negative beta parameter
       c_p     - Positive coefficient
       c_n     - Negative coefficient

       Returns:
       psi     - The computed characteristic function value
       """
    # Computation of Lambda
    Gamma = gamma(-alpha)
    lambda_p = c_p * Gamma * beta_p ** alpha
    lambda_n = c_n * Gamma * beta_n ** alpha
    lambda_tot = lambda_p + lambda_n

    # Define integrands for real and imaginary parts
    def integrand_p_real(z, ui):
        return np.real((z - 1j * ui) ** alpha / ((z ** (alpha + 1))*(b*dt)))

    def integrand_p_imag(z, ui):
        return np.imag((z - 1j * ui) ** alpha / ((z ** (alpha + 1))*(b*dt)))

    def integrand_n_real(z, ui):
        return np.real((z + 1j * ui) ** alpha / ((z ** (alpha + 1))*(b*dt)))

    def integrand_n_imag(z, ui):
        return np.imag((z + 1j * ui) ** alpha / ((z ** (alpha + 1))*(b*dt)))

    # Calculate integrals for each value of u
    int_p_real = np.array(
        [quad(integrand_p_real, beta_p, beta_p * np.exp(b * dt), args = (ui,), limit = 100)[0] for ui in u])
    int_p_imag = np.array(
        [quad(integrand_p_imag, beta_p, beta_p * np.exp(b * dt), args = (ui,), limit = 100)[0] for ui in u])
    int_n_real = np.array(
        [quad(integrand_n_real, beta_n, beta_n * np.exp(b * dt), args = (ui,), limit = 100)[0] for ui in u])
    int_n_imag = np.array(
        [quad(integrand_n_imag, beta_n, beta_n * np.exp(b * dt), args = (ui,), limit = 100)[0] for ui in u])

    phi_J_p = int_p_real + 1j * int_p_imag
    phi_J_n = int_n_real + 1j * int_n_imag

    psi = np.log((np.exp((lambda_tot * dt) * (1 / lambda_tot * (lambda_p * phi_J_p + lambda_n * phi_J_n))) - 1) / (
                np.exp(lambda_tot * dt) - 1))

    return psi
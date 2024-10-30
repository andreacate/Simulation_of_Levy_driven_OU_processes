import numpy as np
from scipy.special import gamma
from scipy.integrate import quad
from scipy.fft import fft
from scipy.stats import poisson, uniform, gamma as gamma_dist, moment
from scipy.interpolate import interp1d
import time
import math
import Cumulants
import CharacteristicFunctions




def integral_via_fft_CDF(M_fft, du, dt, alpha, a, b, beta_p, beta_n, c_p, c_n, gamma_c, scale, flag):
    """
    Function which computes the integral in the Lewis formula for the price
    of a call option via the Fast Fourier Transform

    Parameters:
    M_fft (int):         Power of 2 to compute the number of intervals
    du (float):          Step of the u grid
    dt (float):          Time increment
    alpha (float):       Shape parameter
    a (float):           Parameter a
    b (float):           Mean reversion parameter
    beta_p (float):      Positive beta parameter
    beta_n (float):      Negative beta parameter
    c_p (float):         Positive coefficient
    c_n (float):         Negative coefficient
    gamma_c (float):     Constant shift parameter
    scale (float):       Scaling factor
    flag (int):          Determines which characteristic function to use

    Returns:
    tuple: (I, x) where
           I (numpy array): Value of the integral
           x (numpy array): Moneyness grid
    """

    # FFT parameters:
    N = 2**M_fft
    dx = 2 * np.pi / (du * N)
    x1 = -(N - 1) * dx / 2
    x = np.linspace(x1, -x1, N)
    u1 = -du * (N - 1) / 2
    u = np.linspace(u1, -u1, N)


    if flag == 1:
        integrand = np.exp(CharacteristicFunctions.psi_V(scale * (u + 1j * a), dt, alpha, b, beta_p, beta_n, c_p, c_n)) / (1j * (u + 1j * a))
    elif flag == 2:
        psi_Z = (CharacteristicFunctions.psi_X_TS_OU(scale * (u + 1j * a), alpha, beta_p, beta_n, c_p, c_n, gamma_c) -
                 CharacteristicFunctions.psi_X_TS_OU(scale * (u + 1j * a) * np.exp(-b * dt), alpha, beta_p, beta_n, c_p, c_n, gamma_c))
        integrand = np.exp(psi_Z) / (1j * (u + 1j * a))
    elif flag == 3:
        integrand = np.exp(CharacteristicFunctions.psi_Z_OU_TS(scale * (u + 1j * a), dt, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c)) / (1j * (u + 1j * a))
    else:
        raise ValueError("Invalid flag value")




    # Compute the Fast Fourier Transform
    FFT = fft(integrand * np.exp(-1j * x1 * du * np.arange(N)))

    # Compute the integral by multiplying for the prefactor

    I = np.real(du * np.exp(-1j * u1 * x) * FFT)

    return I, x

def optimal_du_formula(alpha, a, b, c_p, c_n, dt, M, scale, flag):
    """
    omputes the optimal du value based on given parameters.

    Parameters:
        alpha: Stability parameter
        b: Mean reverting parameter
        c_p: Positive c
        c_n: Negative c
        dt:  Time step
        M: Number of time steps
        scale: scaling factor
        flag: model flag

    Returns:
        du: optimal du value
    """
    # Quantities of interest
    N = 2 ** M

    # Compute l based on the flag
    if flag == 2:
        l = -(c_p + c_n) * gamma(-alpha) * np.cos(alpha * np.pi / 2) * (1 - np.exp(-alpha * b * dt)) * scale**alpha
    elif flag == 3:
        l = -(c_p + c_n) * gamma(-alpha) * np.cos(alpha * np.pi / 2) * (1 - np.exp(-alpha * b * dt)) / (alpha * b) * scale**alpha

    # Computation of du
    du = ((2 * np.pi * abs(a)) / (l * N**alpha))**(1 / (alpha + 1))

    return du

def selectCDF(x, discreteCDF, toll):
    """
    Selects a valid segment from the discrete cumulative distribution function (CDF) based on a tolerance.

    Parameters:
    x (array-like): The x-values corresponding to the discrete CDF.
    discreteCDF (array-like): The discrete CDF values.
    toll (float): The tolerance for the CDF condition.

    Returns:
    x_CDF: containing the valid segment of x.
    values_CDF: containing the valid segment of discrete CDF.
    """

    # Check condition
    # Ensure that the difference between consecutive CDF values is within tolerance,
    # and that the CDF values are between 0 and 1.
    check = ((discreteCDF[1:] > (discreteCDF[:-1] - toll)) &
             (discreteCDF[1:] >= 0) &
             (discreteCDF[1:] <= 1))

    # Extend check to include boundaries
    # Add a 0 at the start and end of the check array to consider boundaries
    check = np.concatenate(([0], check, [0]))

    # Find indices where check is zero
    # This will indicate the boundaries of invalid segments
    indexes = np.where(check == 0)[0]

    # Find the segment with the maximum length where check is zero
    # Calculate the differences between consecutive indices
    index_diff = indexes[1:] - indexes[:-1]

    # Find the index of the maximum difference
    indexI = np.argmax(index_diff)

    # Select the range in x and discreteCDF
    # Extract the valid segment using the indices found
    x_CDF = x[indexes[indexI]+1 : indexes[indexI+1]]
    values_CDF = discreteCDF[indexes[indexI]+1 : indexes[indexI+1]]

    return x_CDF, values_CDF


def sampling_from_discrete_CDF(discrete_x, discrete_CDF, N_sample):
    """
    Generates samples from a given discrete cumulative distribution function (CDF) using antithetic variables.

    Parameters:
    discrete_x (array-like): The x-values corresponding to the discrete CDF.
    discrete_CDF (array-like): The discrete CDF values.
    N_sample (int): The number of samples to generate.

    Returns:
    np.ndarray: An array of generated samples.
    """

    # Negative exponential tail
    # Compute parameters for the negative exponential tail fit
    a_n = discrete_CDF[0]
    b_n = np.log(discrete_CDF[1] / discrete_CDF[0]) / (discrete_x[1] - discrete_x[0])

    # Positive exponential tail
    # Compute parameters for the positive exponential tail fit
    a_p = discrete_CDF[-1]
    b_p = np.log((1 - discrete_CDF[-2]) / (1 - discrete_CDF[-1])) / (discrete_x[-2] - discrete_x[-1])

    # Uniform sampling
    # Generate uniform random samples
    u = np.random.rand(N_sample)

    # Ensure the CDF values and corresponding x values are unique
    discrete_CDF_unique, ia = np.unique(discrete_CDF, return_index=True)
    discrete_x_unique = discrete_x[ia]

    def inv_CDF(z):
        """
        Inverse CDF function for interpolation and extrapolation.

        Parameters:
        z (float or array-like): Random values to map through the inverse CDF.

        Returns:
        np.ndarray: The x-values corresponding to the CDF values z.
        """
        # Interpolation within the range of the CDF
        interp_values = interp1d(discrete_CDF_unique, discrete_x_unique, kind='linear', fill_value='extrapolate')(z)

        # Handle values outside the range using exponential tails
        below_range = (np.log(z / a_n) / b_n + discrete_x[0]) * (z <= discrete_CDF[0])
        above_range = (np.log(-(z - 1) / a_p) / b_p + discrete_x[-1]) * (z >= discrete_CDF[-1])

        # Return interpolated values, using tails for out-of-range values
        return np.where(z <= discrete_CDF[0], below_range,
                        np.where(z >= discrete_CDF[-1], above_range, interp_values))

    # Antithetic variables
    # Generate samples using the inverse CDF
    sample_values = np.concatenate((inv_CDF(u), inv_CDF(1 - u)))

    return sample_values




def sim_OU_TS_FinVar_FGMC(x0, b, alpha, beta_p, beta_n, c_p, c_n, gamma_c, T, Nsim, M, M_fft, scale):
    """
    Simulate the OU-TS process with finite variation using FGMC method.

    Parameters:
        x0: Initial condition
        b: Mean reverting parameter
        alpha: Stability parameter
        beta_p: Positive beta
        beta_n: Negative beta
        c_p: Positive c
        c_n: Negative c
        gamma_c: Drift
        T: Total time
        Nsim: Number of simulations
        M: Number of time steps
        M_fft: FFT parameter
        scale: Scaling factor

    Returns:
        X: Simulated paths
        theorCumulantsT: Theoretical cumulants at time T
        simCumulantsT: Simulated cumulants at time T
        theorCumulants_dt: Theoretical cumulants at first time step
        simCumulants_dt: Simulated cumulants at first time step
        logFwd: Log forward in the risk neutral measure
    """

    # Quantities of interest
    a = -0.5 * max(beta_p, beta_n) * ((beta_p > beta_n) - (beta_n > beta_p))
    R = (a < 0)
    dt = T / M
    flagOU_TS_FV = 3

    # Search of optimal du with formula
    du_optimal = optimal_du_formula(alpha, a, b, c_p, c_n, dt, M_fft, scale, flagOU_TS_FV)

    # CDF computation with FFT
    integral_FFT_u, x_grid = integral_via_fft_CDF(M_fft, du_optimal, dt, alpha,
                                                a, b, beta_p, beta_n, c_p, c_n, gamma_c, scale, flagOU_TS_FV)

    # Computation of CDF
    discreteCDF_FFT = R - np.exp(a * x_grid) / np.pi * 0.5 * integral_FFT_u

    # Selecting the discrete CDF
    toll = 1e-7
    x_CDF, values_CDF = selectCDF(x_grid, discreteCDF_FFT, toll)

    # Sampling from the CDF
    X = np.zeros((Nsim, M + 1))
    X[:, 0] = x0
    drift = np.zeros((1, M))

    for ii in range(M):
        sampleValues = sampling_from_discrete_CDF(x_CDF, values_CDF, Nsim // 2)
        X[:, ii + 1] = X[:, ii] * np.exp(-b * dt) + sampleValues / scale
        drift[:,ii] = np.real(-CharacteristicFunctions.psi_Z_OU_TS(-1j, dt * (ii + 1), alpha, b, beta_p, beta_n, c_p, c_n, gamma_c))

    # Computing the log forward in the risk neutral measure
    drift = np.concatenate((np.array([[0]]), drift), axis=1)
    logFwd = X +  np.tile(drift, (Nsim, 1))

    # Computing cumulants
    simCumulantsT = Cumulants.compCumulants(X[:, M]) * 1000
    theorCumulantsT = Cumulants.bctsCumulants(x0, alpha, beta_p, beta_n, c_p, c_n, gamma_c, T, b, flagOU_TS_FV) * 1000

    simCumulants_dt = Cumulants.compCumulants(X[:, 1]) * 1000
    theorCumulants_dt = Cumulants.bctsCumulants(x0, alpha, beta_p, beta_n, c_p, c_n, gamma_c, dt, b, flagOU_TS_FV) * 1000

    return X, theorCumulantsT, simCumulantsT, theorCumulants_dt, simCumulants_dt, logFwd


def sim_TS_OU_FinVar_FGMC(x0, b, alpha, beta_p, beta_n, c_p, c_n, gamma_c, T, Nsim, M, M_fft, scale):
    """
     Simulate the TS-OU process with finite variation using FGMC method.

     Parameters:
         x0: Initial condition
         b: Mean reverting parameter
         alpha: Stability parameter
         beta_p: Positive beta
         beta_n: Negative beta
         c_p: Positive c
         c_n: Negative c
         gamma_c: Drift
         T: Total time
         Nsim: Number of simulations
         M: Number of time steps
         M_fft: FFT parameter
         scale: Scaling factor

     Returns:
         X: Simulated paths
         theorCumulantsT: Theoretical cumulants at time T
         simCumulantsT: Simulated cumulants at time T
         theorCumulants_dt: Theoretical cumulants at first time step
         simCumulants_dt: Simulated cumulants at first time step
         logFwd: Log forward in the risk neutral measure
     """
    ## Quantities of interest
    a = -0.5 * max(beta_p, beta_n) * ((beta_p > beta_n) - (beta_n > beta_p))
    R = a < 0
    dt = T / M
    flagTS_OU = 2

    ## Search of optimal du with formula
    du_optimal = optimal_du_formula(alpha, a, b, c_p, c_n, dt, M_fft, scale, flagTS_OU)

    ## CDF computation with FFT
    integral_FFT_u, x_grid = integral_via_fft_CDF(M_fft, du_optimal, dt, alpha, a, b, beta_p, beta_n, c_p, c_n, gamma_c, scale, flagTS_OU)

    ## Computation of CDF
    discreteCDF_FFT = R - np.exp(a * x_grid) / np.pi * 0.5 * integral_FFT_u

    ## Selecting the discrete CDF
    toll = 1e-9
    x_CDF, values_CDF = selectCDF(x_grid, discreteCDF_FFT, toll)

    ## Sampling from the CDF
    X = np.zeros((Nsim, M + 1))
    X[:, 0] = x0
    drift = np.zeros((1, M))

    for ii in range(M):
        sampleValues = sampling_from_discrete_CDF(x_CDF, values_CDF, Nsim // 2)
        X[:, ii + 1] = X[:, ii] * np.exp(-b * dt) + sampleValues
        drift[:, ii] = np.real(-(CharacteristicFunctions.psi_X_TS_OU(-1j, alpha, beta_p, beta_n, c_p, c_n, gamma_c)
              - CharacteristicFunctions.psi_X_TS_OU(-1j * np.exp(-b * dt*(ii+1)), alpha, beta_p, beta_n, c_p, c_n, gamma_c)))

    ## Computing the log forward in the risk neutral measure
    drift = np.concatenate((np.array([[0]]), drift), axis = 1)
    logFwd = X + np.tile(drift, (Nsim, 1))

    ## Cumulants computation
    # Checking the last time step
    simCumulantsT = Cumulants.compCumulants(X[:, M]) * 1000
    theorCumulantsT = Cumulants.bctsCumulants(x0, alpha, beta_p, beta_n, c_p, c_n, gamma_c, T, b, flagTS_OU) * 1000

    # Checking the first time step
    simCumulants_dt = Cumulants.compCumulants(X[:, 1]) * 1000
    theorCumulants_dt = Cumulants.bctsCumulants(x0, alpha, beta_p, beta_n, c_p, c_n, gamma_c, dt, b, flagTS_OU) * 1000

    return X, theorCumulantsT, simCumulantsT, theorCumulants_dt, simCumulants_dt, logFwd

def sim_OU_TS_FinAct_FGMC(x0, b, alpha, beta_p, beta_n, c_p, c_n, gamma_c, T, Nsim, M, M_fft, scale):
    """
    Simulate the OU-TS process with finite activity using FGMC method.

    Parameters:
        x0: Initial condition
        b: Mean reverting parameter
        alpha: Stability parameter
        beta_p: Positive beta
        beta_n: Negative beta
        c_p: Positive c
        c_n: Negative c
        gamma_c: Drift
        T: Total time
        Nsim: Number of simulations
        M: Number of time steps
        M_fft: FFT parameter
        scale: Scaling factor

    Returns:
        X: Simulated paths
        theorCumulantsT: Theoretical cumulants at time T
        simCumulantsT: Simulated cumulants at time T
        theorCumulants_dt: Theoretical cumulants at first time step
        simCumulants_dt: Simulated cumulants at first time step
        logFwd: Log forward in the risk neutral measure
    """

    # Quantities of interest
    a = -0.25 * max(beta_p, beta_n) * ((beta_p > beta_n) - (beta_n > beta_p))
    R = a < 0
    dt = T / M
    flagOU_TS_FA = 1

    # Lambda computation
    Gamma = gamma(-alpha)
    lambda_p = c_p * Gamma * beta_p**alpha
    lambda_n = c_n * Gamma * beta_n**alpha
    lambda_tot = lambda_p + lambda_n

    # Sampling of the Bernoulli
    U = np.random.rand(Nsim, M)
    Bern = U < (1 - np.exp(-lambda_tot * dt))

    # CDF computation with FFT
    du_optimal = 0.4
    integral_FFT_u, x_grid = integral_via_fft_CDF(M_fft, du_optimal, dt, alpha,
                                                a, b, beta_p, beta_n, c_p, c_n, gamma_c, scale, flagOU_TS_FA)

    # Computation of CDF
    discreteCDF_FFT = R - np.exp(a * x_grid) / np.pi * 0.5 * integral_FFT_u

    # Selecting the discrete CDF
    toll = 0
    x_CDF, values_CDF = selectCDF(x_grid, discreteCDF_FFT, toll)

    # Computation of X
    X = np.zeros((Nsim, M + 1))
    X[:, 0] = x0
    drift = np.zeros((1, M))

    for ii in range(M):
        sampleValues = sampling_from_discrete_CDF(x_CDF, values_CDF, Nsim // 2)
        X[:, ii + 1] = X[:, ii] * np.exp(-b * dt) + sampleValues * Bern[:, ii] / scale
        drift[:, ii] = np.real(
            -CharacteristicFunctions.psi_Z_OU_TS(-1j, dt * (ii + 1), alpha, b, beta_p, beta_n, c_p, c_n, gamma_c))

    # Computing the log forward in the risk neutral measure
    drift = np.concatenate((np.array([[0]]), drift), axis = 1)
    logFwd = X + np.tile(drift, (Nsim, 1))

    # Cumulants computation
    simCumulantsT = Cumulants.compCumulants(X[:, M]) * 1000
    theorCumulantsT = Cumulants.bctsCumulants(x0, alpha, beta_p, beta_n, c_p, c_n, gamma_c, T, b, flagOU_TS_FA) * 1000

    simCumulants_dt = Cumulants.compCumulants(X[:, 1]) * 1000
    theorCumulants_dt = Cumulants.bctsCumulants(x0, alpha, beta_p, beta_n, c_p, c_n, gamma_c, dt, b, flagOU_TS_FA) * 1000

    return X, theorCumulantsT, simCumulantsT, theorCumulants_dt, simCumulants_dt, logFwd





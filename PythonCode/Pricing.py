import numpy as np
from scipy.interpolate import interp1d
from scipy.fft import fft
from scipy.stats import norm
import CharacteristicFunctions
import ExactDecomposition
import FastGeneralMonteCarlo
import numpy as np



def integral_via_fft(f, M, x1):
    """
    Function which computes the integral in the Lewis formula for the price
    of a call option via the Fast Fourier Transform.

    Parameters:
    f: function to be Fourier transformed
    M: power of 2 to compute the number of intervals
    x1: parameter related to the integration range

    Returns:
    I: value of the integral
    x: moneyness grid
    """
    # FFT parameters
    N = 2 ** M
    dx = -2 * x1 / (N - 1)
    x = np.linspace(x1, -x1, N)
    du = 2 * np.pi / (dx * N)
    u1 = -du * (N - 1) / 2
    u = np.linspace(u1, -u1, N)

    # Compute the Fast Fourier Transform
    FFT = fft(f(u) * np.exp(-1j * x1 * du * np.arange(N)))

    # Compute the integral by multiplying for the prefactor
    I = np.real(du * np.exp(-1j * u1 * x) * FFT)

    return I, x


def price_european_lewis_fft(S0, b, alpha, beta_p, beta_n, c_p, c_n, gamma_c, T, moneyness, r, flag, M_fft):
    """
    Computation of European Call Price via Lewis formula

    INPUT:
    S0:        initial underlying condition
    b:         mean reverting parameter
    alpha:     stability parameter
    beta_p:    positive beta
    beta_n:    negative beta
    c_p:       positive c
    c_n:       negative c
    gamma_c:   drift
    T:         time to maturity
    moneyness: vector of moneyness
    r:         rate
    flag:      1 -> OU-CTS Finite Activity
              2 -> CTS-OU Finite Variation
               3 -> OU-CTS Finite Variation
    M_fft:     parameter of FFT

    USES:
    function psi_Z_OU_TS(u, dt, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c)
    function psi_X_TS_OU(u, alpha, beta_p, beta_n, c_p, c_n, gamma_c)
    function integralViaFFT(f, M, x1)
    """

    if flag == 1:
        drift = -CharacteristicFunctions.psi_Z_OU_TS(-1j, T, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c)
        integrand = lambda u: np.exp(CharacteristicFunctions.psi_Z_OU_TS(-u - 0.5j, T, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c)) / (
                    u ** 2 + 0.25) * np.exp(1j * (-u - 0.5j) * drift)

    elif flag == 2:
        drift = -(CharacteristicFunctions.psi_X_TS_OU(-1j, alpha, beta_p, beta_n, c_p, c_n, gamma_c)
                  - CharacteristicFunctions.psi_X_TS_OU(-1j * np.exp(-b * T), alpha,beta_p, beta_n, c_p, c_n,gamma_c))
        integrand = lambda u: (np.exp(CharacteristicFunctions.psi_X_TS_OU(-u - 0.5j, alpha, beta_p, beta_n, c_p, c_n, gamma_c) -
                                      CharacteristicFunctions.psi_X_TS_OU((-u - 0.5j) * np.exp(-b * T), alpha, beta_p, beta_n, c_p, c_n,
                                                  gamma_c)) /
                               (u ** 2 + 0.25) * np.exp(1j * (-u - 0.5j) * drift))

    elif flag == 3:
        drift = -CharacteristicFunctions.psi_Z_OU_TS(-1j, T, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c)
        integrand = lambda u: np.exp(CharacteristicFunctions.psi_Z_OU_TS(-u - 0.5j, T, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c)) / (
                    u ** 2 + 0.25) * np.exp(1j * (-u - 0.5j) * drift)

    integral_fft, x_grid = integral_via_fft(integrand, M_fft, moneyness[0] * 1000)
    integral_lewis = interp1d(x_grid, integral_fft, kind = 'cubic')(moneyness)

    # European price
    fwd = S0 * np.exp(r * T)
    prices = np.exp(-r * T) * fwd * (1 - np.exp(-moneyness / 2) * np.real(integral_lewis) / (2 * np.pi))

    return prices



def priceEuropean(x0, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c, T, M, M_fft, Nsim, r, S0, moneyness, scale, model, method):
    """
    Computation of European Call Price via Lewis formula.

    Parameters:
    S0: Initial stock price
    b: Parameter for the model
    alpha: Parameter for the model
    beta_p: Parameter for the model
    beta_n: Parameter for the model
    c_p: Parameter for the model
    c_n: Parameter for the model
    gamma_c: Parameter for the model
    T: Time to maturity
    moneyness: Array of moneyness values
    r: Risk-free rate
    flag: Model selection flag
    M_fft: Power of 2 for FFT intervals

    Returns:
    prices: Computed European call prices
    """
    dt = T / M

    # Simulation
    if model == 1:
        if method == 1:
            _, _, _, _, _, logFwd = ExactDecomposition.sim_OU_TS_FinAct_ED(x0, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c, T, M, Nsim)
        else:
            _, _, _, _, _, logFwd = FastGeneralMonteCarlo.sim_OU_TS_FinAct_FGMC(x0, b, alpha, beta_p, beta_n, c_p, c_n, gamma_c, T, Nsim, M, M_fft, scale)
    elif model == 2:
        if method == 1:
            _, _, _, _, _, logFwd = ExactDecomposition.sim_TS_OU_FinVar_ED(x0, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c, T, M, Nsim)
        else:
            _, _, _, _, _, logFwd = FastGeneralMonteCarlo.sim_TS_OU_FinVar_FGMC(x0, b, alpha, beta_p, beta_n, c_p, c_n, gamma_c, T, Nsim, M, M_fft, scale)
    elif model == 3:
        if method == 1:
            _, _, _, _, _, logFwd = ExactDecomposition.sim_OU_TS_FinVar_ED(x0, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c, T, M, Nsim)
        else:
            _, _, _, _, _, logFwd = FastGeneralMonteCarlo.sim_OU_TS_FinVar_FGMC(x0, b, alpha, beta_p, beta_n, c_p, c_n, gamma_c, T, Nsim, M, M_fft, scale)

    # Ensure logFwd contains no complex values
    logFwd = np.real(logFwd)

    # Forward
    fwd = S0 * np.exp(r * T)

    # Underlying
    underlying_sim = fwd * np.exp(logFwd[:, -1])

    # Strikes
    K = fwd * np.exp(-moneyness)

    # Discounted payoff
    payoff = np.maximum(np.tile(underlying_sim, (len(K), 1)).T - K, 0)
    disc_payoff = payoff * np.exp(-r * T)

    # Prices
    prices, std_dev = np.mean(disc_payoff, axis=0), np.std(disc_payoff, axis=0)
    conf_intervals = norm.interval(0.95, loc=prices, scale=std_dev / np.sqrt(Nsim))

    return prices, std_dev, conf_intervals






def priceAmericanLS(S, T, K, r):
    """
    Prices American options using Least-Squares Monte Carlo (LSMC) method.

    Parameters:
    S (numpy.ndarray): Simulated stock prices (N x M+1) where N is the number of paths and M is the number of time steps.
    T (float): Time to maturity.
    K (numpy.ndarray): Strike prices.
    r (float): Risk-free interest rate.

    Returns:
    numpy.ndarray: Prices of American options for each strike price in K.
    """
    N, M = S.shape
    M -= 1  # Adjust M to be the number of time steps
    dt = T / M
    exerciseTime = M * np.ones(N, dtype=int)
    prices = np.zeros(len(K))
    conf_intervals = np.zeros((2, len(K)))

    for jj in range(len(prices)):
        payoffAM = np.maximum(0, S[:, -1] - K[jj])

        for ii in range(M, 0, -1):
            IntrValue = np.maximum(0, S[:, ii] - K[jj])
            indexITM = np.where(IntrValue > 0)[0]

            if len(indexITM) == 0:
                continue

            b = payoffAM[indexITM] * np.exp(-r * dt * (exerciseTime[indexITM] - ii + 1))
            A = np.vstack([np.ones(len(indexITM)), S[indexITM, ii], S[indexITM, ii] ** 2]).T

            weights = np.linalg.lstsq(A, b, rcond=None)[0]

            CV = A @ weights
            IV = IntrValue[indexITM]

            indexEarlyEx = np.where(IV > CV)[0]

            payoffAM[indexITM[indexEarlyEx]] = IV[indexEarlyEx]
            exerciseTime[indexITM[indexEarlyEx]] = ii - 1

        discPrice = payoffAM * np.exp(-r * dt * exerciseTime)
        # Prices
        mean_price, std_dev = np.mean(discPrice), np.std(discPrice)
        # Calculate the confidence interval using norm.ppf
        z_score = norm.ppf(0.975)  # for a 95% confidence interval (two-tailed)
        margin_of_error = z_score * std_dev / np.sqrt(N)

        prices[jj] = mean_price
        conf_intervals[0, jj] = mean_price - margin_of_error  # Lower bound
        conf_intervals[1, jj] = mean_price + margin_of_error  # Upper bound

    return prices, conf_intervals


def priceAmerican(x0, S, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c, T, M, M_fft, Nsim, K, r, scale, model, method):
    """
    Computation of American Call Option Price

    Parameters:
    x0: Initial underlying asset price
    alpha: Stability parameter
    b: Mean reverting parameter
    beta_p: Positive beta
    beta_n: Negative beta
    c_p: Positive c
    c_n: Negative c
    gamma_c: Drift term
    T: Time to maturity
    M: Number of time steps for simulation
    M_fft: Parameter for FFT (if applicable)
    Nsim: Number of simulations
    K: Strike prices (array)
    r: Risk-free interest rate
    scale: Scaling factor (if applicable for FGMC method)
    model: Model selection flag (1: OU-CTS Finite Activity, 2: CTS-OU Finite Variation, 3: OU-CTS Finite Variation)
    method: Method selection flag (1: Exact Decomposition, 2: FGMC)

    Returns:
    price: American call option price
    conf_int: Confidence intervals for the option price (optional)

    Uses:
    * sim_OU_TS_FinAct_ED(x0, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c, T, M, Nsim) (for OU-CTS Finite Activity with Exact Decomposition)
    * sim_OU_TS_FinAct_FGMC(x0, b, alpha, beta_p, beta_n, c_p, c_n, gamma_c, T, Nsim, M, M_fft, scale) (for OU-CTS Finite Activity with FGMC)
    * (similar function definitions for other models and methods)
    * priceAmericanLS(S, T, K, r) (potentially used for Longstaff-Schwartz calculation)
    """
    dt = T / M  # Time step size

    # Simulation
    if model == 1:
        # OU-TS Finite Activity model
        if method == 1:
            # Exact decomposition method
            _, _, _, _, _, logFwd = ExactDecomposition.sim_OU_TS_FinAct_ED(x0, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c, T, M, Nsim)
        else:
            # Fast General Monte Carlo method
            _, _, _, _, _, logFwd = FastGeneralMonteCarlo.sim_OU_TS_FinAct_FGMC(x0, b, alpha, beta_p, beta_n, c_p, c_n, gamma_c, T, Nsim, M,
                                                          M_fft, scale)
    elif model == 2:
        # TS-OU Finite Variation model
        if method == 1:
            # Exact decomposition method
            _, _, _, _, _, logFwd = ExactDecomposition.sim_TS_OU_FinVar_ED(x0, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c, T, M, Nsim)
        else:
            # Fast General Monte Carlo method
            _, _, _, _, _, logFwd = FastGeneralMonteCarlo.sim_TS_OU_FinVar_FGMC(x0, b, alpha, beta_p, beta_n, c_p, c_n, gamma_c, T, Nsim, M,
                                                          M_fft, scale)
    elif model == 3:
        # OU-TS Finite Variation model
        if method == 1:
            # Exact decomposition method
            _, _, _, _, _, logFwd = ExactDecomposition.sim_OU_TS_FinVar_ED(x0, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c, T, M, Nsim)
        else:
            # Fast General Monte Carlo method
            _, _, _, _, _, logFwd = FastGeneralMonteCarlo.sim_OU_TS_FinVar_FGMC(x0, b, alpha, beta_p, beta_n, c_p, c_n, gamma_c, T, Nsim, M,
                                                          M_fft, scale)

    # Price computation
    fwd = S * np.exp(r*T)

    F = fwd * np.exp(logFwd)  # Forward prices

    S = F * np.exp(-r * dt * np.arange(M, -1, -1))  # Discounted stock prices

    price, conf_intervals = priceAmericanLS(S, T, K, r)

    return price, conf_intervals


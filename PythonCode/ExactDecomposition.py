
import numpy as np
from scipy.special import gamma
from scipy.stats import poisson, uniform, gamma as gamma_dist, moment
import math
import time
import matplotlib.pyplot as plt
import Cumulants
import CharacteristicFunctions





def sim_OU_TS_FinAct_ED(x0, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c, T, M, Nsim):
    """
    Simulation of OU-TS Finite Activity with Exact Decomposition following
    the algorithm 1 in Sabino [3]

    INPUT
    x0:       initial condition
    alpha:    stability parameter
    b:        mean reverting parameter
    beta_p:   positive beta
    beta_n:   negative beta
    c_p:      positive c
    c_n:      negative c
    gamma_c:  drift
    T:        time to maturity
    M:        number of steps
    Nsim:     number of simulations

    OUTPUT
    X:                   matrix with simulations
    cumulants_process:   theoretical cumulants of the process
    cumulants_sim:       cumulants of the simulations
    logFwd:              log of fwd prices

    USES
    function theorCumulants(x0, alpha, beta, c, dt, b, flag)
    function bctsCumulants(x0, alpha, beta_p, beta_n, c_p, c_n, gamma_c, dt, b, flag)
    function compCumulants(vec)
    function psi_Z_OU_TS(u, dt, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c)
    """
    flag_OU_TS_FA = 1
    X = np.zeros((Nsim, M+1))
    dt = T / M
    drift = np.zeros((1, M))

    # Lambda computation
    Gamma = gamma(-alpha)
    a = np.exp(-b * dt)
    lambda_p = c_p * Gamma * beta_p**alpha
    lambda_n = c_n * Gamma * beta_n**alpha

    # Simulations
    poisson_p = poisson.rvs(lambda_p * dt, size=(Nsim, M))
    poisson_n = poisson.rvs(lambda_n * dt, size=(Nsim, M))

    for ii in range(M):
        max_p = np.max(poisson_p[:, ii])
        max_n = np.max(poisson_n[:, ii])

        uniform_p = uniform.rvs(size=(Nsim, max_p))
        uniform_n = uniform.rvs(size=(Nsim, max_n))

        beta_hat_p = beta_p * np.exp(b * uniform_p * dt)
        beta_hat_n = beta_n * np.exp(b * uniform_n * dt)

        j_hat_p = gamma_dist.rvs(-alpha, scale=1./beta_hat_p)
        j_hat_n = gamma_dist.rvs(-alpha, scale=1./beta_hat_n)

        count_matrix_p = np.tile(np.arange(1, max_p + 1), (Nsim, 1))
        count_matrix_n = np.tile(np.arange(1, max_n + 1), (Nsim, 1))

        rep_matrix_p = np.tile(poisson_p[:, ii], (max_p, 1)).T
        rep_matrix_n = np.tile(poisson_n[:, ii], (max_n, 1)).T

        find_p = rep_matrix_p >= count_matrix_p
        find_n = rep_matrix_n >= count_matrix_n

        X2_p = np.sum(j_hat_p * find_p, axis=1)
        X2_n = np.sum(j_hat_n * find_n, axis=1)

        cumulants_p = Cumulants.theorCumulants(x0, alpha, beta_p, c_p, dt, b, flag_OU_TS_FA)
        cumulants_n = Cumulants.theorCumulants(x0, alpha, beta_n, c_n, dt, b, flag_OU_TS_FA)

        X[:, ii+1] = gamma_c * dt + X[:, ii] * a + X2_p - cumulants_p[0] - (X2_n - cumulants_n[0])
        drift[:, ii] = np.real(
            -CharacteristicFunctions.psi_Z_OU_TS(-1j, dt * (ii + 1), alpha, b, beta_p, beta_n, c_p, c_n, gamma_c))

    ## Cumulants computation
    # Checking the last time step
    simCumulantsT = Cumulants.compCumulants(X[:, M]) * 1000
    theorCumulantsT = Cumulants.bctsCumulants(x0, alpha, beta_p, beta_n, c_p, c_n, gamma_c, T, b, flag_OU_TS_FA) * 1000

    # Checking the first time step
    simCumulants_dt = Cumulants.compCumulants(X[:, 1]) * 1000
    theorCumulants_dt = Cumulants.bctsCumulants(x0, alpha, beta_p, beta_n, c_p, c_n, gamma_c, dt, b, flag_OU_TS_FA) * 1000

    # Computing the log forward in the risk-neutral measure
    drift = np.concatenate((np.array([[0]]), drift), axis = 1)
    logFwd = X + np.tile(drift, (Nsim, 1))

    return X, theorCumulantsT, simCumulantsT, theorCumulants_dt, simCumulants_dt, logFwd




def simulateTS(alpha, beta, c):
    # Acceptance Rejection Simulation
    accettato = False

    while not accettato:
        # [1] generate a stable random variable via Zolotarev
        Us = (np.random.rand() - 0.5) * np.pi
        Es = -np.log(np.random.rand())  # generating exponential r.v. (exprnd is slow!)
        S = (-c * gamma(-alpha))**(1/alpha) * \
            np.sin(alpha * Us + 0.5*np.pi*alpha) / np.cos(Us)**(1/alpha) * \
            (np.cos((1-alpha)*Us - 0.5*np.pi*alpha) / Es)**((1-alpha)/alpha)

        # [2] generate a uniform U
        U = np.random.rand()

        # [3] A/R condition
        accettato = (U <= np.exp(-beta*S))

    # Corrections
    # add "drift" part of the characteristic function:
    # the one associated with the third term of the Lévy-Khintchine
    S_mean = c * beta**(alpha-1) * gamma(-alpha+1)

    return S, S_mean

def simulateV(alpha, a, nSim):
    L = 100

    # Define the pdf W, the modification of the original pdf V which is
    # convex and increasing
    def pdf_W(x):
        return -a**alpha * np.log(a**alpha) / \
            (1 - a**alpha + a**alpha * np.log(a**alpha)) * \
            (a**(-alpha*x) - 1)

    # Generate a grid with L+1 ticks
    w_grid = np.linspace(0, 1, L+1)

    # Compute the corresponding values for the PDF
    fw_grid = np.zeros_like(w_grid)
    for i in range(len(w_grid)):
        fw_grid[i] = pdf_W(w_grid[i])

    # Define the new rescaled quantities
    ql_grid = 0.5 * (fw_grid[:-1] + fw_grid[1:]) / L
    GL = np.sum(ql_grid)
    pl_grid = ql_grid / GL
    cumul_pl_grid = np.cumsum(pl_grid)

    # A/R simulation
    v_vector = np.zeros(nSim)

    for i in range(nSim):
        accettato = False

        while not accettato:
            # Draw a rv s, i.e., to identify a small interval I_s
            s = np.where(cumul_pl_grid > np.random.rand())[0][0]

            # Simulate the rv g_s, corresponding to interval I_s
            coeff_ang = (fw_grid[s+1] - fw_grid[s]) / (w_grid[s+1] - w_grid[s])
            y = w_grid[s] - fw_grid[s] / coeff_ang + \
                np.sqrt(fw_grid[s]**2 / coeff_ang**2 + 2 / coeff_ang * ql_grid[s] * np.random.rand())

            # Simulate the uniform rv "u", for the A/R method
            u = np.random.rand()

            # A/R condition
            fw_y = pdf_W(y)
            gl_y = coeff_ang * (y - w_grid[s]) + fw_grid[s]

            # Exit clause
            accettato = (u <= (fw_y / gl_y))

        # Use the inverse transformation to obtain the original rv V
        v_vector[i] = a**(-y)

    return v_vector


def sim_OU_TS_FinVar_ED(x0, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c, T, M, Nsim):
    """
    Simulation of OU-TS Finite Variation with Exact Decomposition following
    the algorithm 2 in Sabino [4]

    INPUT
    x0:       initial condition
    alpha:    stability parameter
    b:        mean reverting parameter
    beta_p:   positive beta
    beta_n:   negative beta
    c_p:      positive c
    c_n:      negative c
    gamma_c:  drift
    T:        time to maturity
    M:        number of steps
    Nsim:     number of simulations

    OUTPUT
    X:                   matrix with simulations
    cumulants_process:   theoretical cumulants of the process
    cumulants_sim:       cumulants of the simulations
    logFwd:              log of fwd prices

    USES
    function simulateTS(alpha, beta, c)
    function simulateV(alpha, a, nSim)
    function theorCumulants(x0, alpha, beta, c, dt, b, flag)
    function bctsCumulants(x0, alpha, beta_p, beta_n, c_p, c_n, gamma_c, dt, b, flag)
    function compCumulants(vec)
    function psi_Z_OU_TS(u, dt, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c)
    """

    # Quantities of interest
    flag_OU_TS_FV = 3
    X = np.zeros((Nsim, M+1))
    dt = T / M
    drift = np.zeros((1, M))

    # Lambda computation
    Gamma = gamma(1 - alpha)
    a = np.exp(-b * dt)
    lambda_p = (c_p * beta_p ** alpha * Gamma / (b * alpha ** 2 * a ** alpha) *
                (1 - a ** alpha + a ** alpha * np.log(a ** alpha)))
    lambda_n = (c_n * beta_n ** alpha * Gamma / (b * alpha ** 2 * a ** alpha) *
                (1 - a ** alpha + a ** alpha * np.log(a ** alpha)))

    # Simulations
    poisson_p = np.random.poisson(lambda_p, size=(Nsim, M))
    poisson_n = np.random.poisson(lambda_n, size=(Nsim, M))

    for ii in range(M):
        X1_p = np.zeros(Nsim)
        X1_n = np.zeros(Nsim)

        for jj in range(Nsim):
            X1_p[jj], _ = simulateTS(alpha, beta_p / a, c_p * (1 - a ** alpha) / (alpha * b))
            X1_n[jj], _ = simulateTS(alpha, beta_n / a, c_n * (1 - a ** alpha) / (alpha * b))

        max_p = np.max(poisson_p[:, ii])
        max_n = np.max(poisson_n[:, ii])

        V_p = np.zeros((Nsim, max_p))
        V_n = np.zeros((Nsim, max_n))

        for kk in range(max_p):
            V_p[:, kk] = simulateV(alpha, a, Nsim)

        for kk in range(max_n):
            V_n[:, kk] = simulateV(alpha, a, Nsim)

        beta_hat_p = beta_p * V_p
        beta_hat_n = beta_n * V_n

        j_hat_p = np.random.gamma((1 - alpha) * np.ones(beta_hat_p.shape), 1. / beta_hat_p)
        j_hat_n = np.random.gamma((1 - alpha) * np.ones(beta_hat_n.shape), 1. / beta_hat_n)

        count_matrix_p = np.tile(np.arange(1, max_p + 1), (Nsim, 1))
        count_matrix_n = np.tile(np.arange(1, max_n + 1), (Nsim, 1))

        rep_matrix_p = np.tile(poisson_p[:, ii].reshape(-1, 1), (1, max_p))
        rep_matrix_n = np.tile(poisson_n[:, ii].reshape(-1, 1), (1, max_n))

        find_p = rep_matrix_p >= count_matrix_p
        find_n = rep_matrix_n >= count_matrix_n

        cumulants_p = Cumulants.theorCumulants(x0, alpha, beta_p, c_p, dt, b, flag_OU_TS_FV)
        cumulants_n = Cumulants.theorCumulants(x0, alpha, beta_n, c_n, dt, b, flag_OU_TS_FV)

        X2_p = np.sum(j_hat_p * find_p, axis=1)
        X2_n = np.sum(j_hat_n * find_n, axis=1)

        X[:, ii + 1] = (gamma_c * dt + a * X[:, ii] + X1_p + X2_p - cumulants_p[0] -
                        (X1_n + X2_n - cumulants_n[0]))
        drift[:, ii] = np.real(
            -CharacteristicFunctions.psi_Z_OU_TS(-1j, dt * (ii + 1), alpha, b, beta_p, beta_n, c_p, c_n, gamma_c))

    ## Cumulants computation
    # Checking the last time step
    simCumulantsT = Cumulants.compCumulants(X[:, M]) * 1000
    theorCumulantsT = Cumulants.bctsCumulants(x0, alpha, beta_p, beta_n, c_p, c_n, gamma_c, T, b, flag_OU_TS_FV) * 1000

    # Checking the first time step
    simCumulants_dt = Cumulants.compCumulants(X[:, 1]) * 1000
    theorCumulants_dt = Cumulants.bctsCumulants(x0, alpha, beta_p, beta_n, c_p, c_n, gamma_c, dt, b, flag_OU_TS_FV) * 1000


    # Computing the log forward in the risk-neutral measure
    drift = np.concatenate((np.array([[0]]), drift), axis = 1)
    logFwd = X + np.tile(drift, (Nsim, 1))

    return X, theorCumulantsT, simCumulantsT, theorCumulants_dt, simCumulants_dt, logFwd


def sim_TS_OU_FinVar_ED(x0, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c, T, M, Nsim):
    """
    Simulation of TS-OU Finite Variation with Exact Decomposition following
    the algorithm 1 in [4] Sabino

    INPUT
    x0:       initial condition
    alpha:    stability parameter
    b:        mean reverting parameter
    beta_p:   positive beta
    beta_n:   negative beta
    c_p:      positive c
    c_n:      negative c
    gamma_c:  drift
    T:        time to maturity
    M:        number of steps
    Nsim:     number of simulations

    OUTPUT
    X:                   matrix with simulations
    cumulants_process:   theoretical cumulants of the process
    cumulants_sim:       cumulants of the simulations
    logFwd:              log of fwd prices

    USES
    function simulateTS(alpha, beta, c)
    function theorCumulants(x0, alpha, beta, c, dt, b, flag)
    function bctsCumulants(x0, alpha, beta_p, beta_n, c_p, c_n, gamma_c, dt, b, flag)
    function compCumulants(vec)
    """

    flag_TS_OU = 2
    X = np.zeros((Nsim, M+1))
    dt = T/M
    drift = np.zeros((1, M))

    Gamma = gamma(1-alpha)
    a = np.exp(-b*dt)
    lambda_p = c_p * Gamma * beta_p**alpha/alpha * (1-a**alpha)
    lambda_n = c_n * Gamma * beta_n**alpha/alpha * (1-a**alpha)

    poisson_p = np.random.poisson(lambda_p, (Nsim, M))
    poisson_n = np.random.poisson(lambda_n, (Nsim, M))

    for ii in range(M):
        X1_p = np.zeros(Nsim)
        X1_n = np.zeros(Nsim)

        for jj in range(Nsim):
            X1_p[jj], _ = simulateTS(alpha, beta_p, c_p*(1-a**alpha))
            X1_n[jj], _ = simulateTS(alpha, beta_n, c_n*(1-a**alpha))

        max_p = np.max(poisson_p[:,ii])
        max_n = np.max(poisson_n[:,ii])

        uniform_p = np.random.rand(Nsim, max_p)
        uniform_n = np.random.rand(Nsim, max_n)

        V_p = (1 + (a**(-alpha)-1) * uniform_p)**(1/alpha)
        V_n = (1 + (a**(-alpha)-1) * uniform_n)**(1/alpha)

        beta_hat_p = beta_p * V_p
        beta_hat_n = beta_n * V_n

        j_hat_p = np.random.gamma((1-alpha)*np.ones_like(uniform_p), 1./beta_hat_p)
        j_hat_n = np.random.gamma((1-alpha)*np.ones_like(uniform_n), 1./beta_hat_n)

        count_matrix_p = np.tile(np.arange(1, max_p+1), (Nsim, 1))
        count_matrix_n = np.tile(np.arange(1, max_n+1), (Nsim, 1))

        rep_matrix_p = np.tile(poisson_p[:,ii], (max_p, 1)).T
        rep_matrix_n = np.tile(poisson_n[:,ii], (max_n, 1)).T

        find_p = rep_matrix_p >= count_matrix_p
        find_n = rep_matrix_n >= count_matrix_n

        cumulants_p = Cumulants.theorCumulants(x0, alpha, beta_p, c_p, dt, b, flag_TS_OU)
        cumulants_n = Cumulants.theorCumulants(x0, alpha, beta_n, c_n, dt, b, flag_TS_OU)

        X2_p = np.sum(j_hat_p*find_p, axis=1)
        X2_n = np.sum(j_hat_n*find_n, axis=1)

        X[:,ii+1] = gamma_c * dt + a * X[:,ii] + X1_p + X2_p - cumulants_p[0] - (X1_n + X2_n - cumulants_n[0])
        drift[:, ii] = np.real(-(CharacteristicFunctions.psi_X_TS_OU(-1j, alpha, beta_p, beta_n, c_p, c_n, gamma_c)
                                 - CharacteristicFunctions.psi_X_TS_OU(-1j * np.exp(-b * dt * (ii + 1)), alpha, beta_p,
                                                                       beta_n, c_p, c_n, gamma_c)))

    ## Cumulants computation
    # Checking the last time step
    simCumulantsT = Cumulants.compCumulants(X[:, M]) * 1000
    theorCumulantsT = Cumulants.bctsCumulants(x0, alpha, beta_p, beta_n, c_p, c_n, gamma_c, T, b, flag_TS_OU) * 1000

    # Checking the first time step
    simCumulants_dt = Cumulants.compCumulants(X[:, 1]) * 1000
    theorCumulants_dt = Cumulants.bctsCumulants(x0, alpha, beta_p, beta_n, c_p, c_n, gamma_c, dt, b, flag_TS_OU) * 1000

    # Computing the log forward in the risk-neutral measure
    drift = np.concatenate((np.array([[0]]), drift), axis = 1)
    logFwd = X + np.tile(drift, (Nsim, 1))

    return X, theorCumulantsT, simCumulantsT, theorCumulants_dt, simCumulants_dt, logFwd


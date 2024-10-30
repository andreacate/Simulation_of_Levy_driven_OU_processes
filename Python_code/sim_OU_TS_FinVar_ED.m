function [X, theorCumulantsT, simCumulantsT, theorCumulants_dt, simCumulants_dt, logFwd] ...
          = sim_OU_TS_FinVar_ED(x0, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c, T, M, Nsim)
% Simulation of OU-TS Finite Variation with Exact Decomposition following
% the algorithm 2 in Sabino [4] 
%
% INPUT
% x0:       initial condition
% alpha:    stability parameter    
% b:        mean reverting parameter
% beta_p:   positive beta
% beta_n:   negative beta
% c_p:      positive c
% c_n:      negative c
% gamma_c:  drift
% T:        time to maturity
% M:        number of steps
% Nsim:     number of simulations
%
% OUTPUT
% X:                   matrix with simulations
% cumulants_process:   theoretical cumulants of the process
% cumulants_sim:       cumulants of the simulations
% logFwd:              log of fwd prices
%
% USES
% function simulateTS(alpha, beta, c)
% function simulateV(alpha, a, nSim)
% function theorCumulants(X0, alpha, beta, c, dt, b, flag)
% function bctsCumulants(X0, alpha, beta_p, beta_n, c_p, c_n, gamma_c, dt, b, flag)
% function compCumulants(vec)
% function psi_Z_OU_TS(u, dt, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c)

    %%  Quantities of interest
    
    flag_OU_TS_FV = 3;  % Flag for OU-TS Finite Variation
    X = zeros(Nsim, M+1);  % Matrix to store simulations
    X(:,1) = x0 * ones(Nsim, 1);  % Set initial condition
    drift = zeros(1, M);  % Initialize drift vector
    dt = T / M;  % Time step size

    % Lambda computation
    Gamma = gamma(1 - alpha);  
    a = exp(-b * dt);  
    lambda_p = c_p * beta_p^alpha * Gamma / (b * alpha^2 * a^alpha) * (1 - a^alpha + a^alpha * log(a^alpha));  
    lambda_n = c_n * beta_n^alpha * Gamma / (b * alpha^2 * a^alpha) * (1 - a^alpha + a^alpha * log(a^alpha));  

    %% Simulations
    
    % Poisson simulation
    poisson_p = random('Poisson', lambda_p, Nsim, M);  % Positive Poisson process
    poisson_n = random('Poisson', lambda_n, Nsim, M);  % Negative Poisson process

    % Process simulation
    for ii = 1:M
        % Inizialization of the vector
        X1_p = zeros(Nsim, 1);  % Positive TS increments
        X1_n = zeros(Nsim, 1);  % Negative TS increments
        
        % Simulating Tempered stable vectors
        [X1_p, ~] = simulateTS(alpha, beta_p/a, c_p * (1 - a^alpha) / (alpha * b), Nsim);  % Simulate positive TS increment
        [X1_n, ~] = simulateTS(alpha, beta_n/a, c_n * (1 - a^alpha) / (alpha * b), Nsim);  % Simulate positive TS increment
        

        % Compound Poisson simulation        
        max_p = max(poisson_p(:, ii));  % Max number of positive jumps
        max_n = max(poisson_n(:, ii));  % Max number of negative jumps

        V_p = zeros(Nsim, max_p);  % Positive jump sizes
        V_n = zeros(Nsim, max_n);  % Negative jump sizes

        for kk = 1:max_p
            V_p(:, kk) = simulateV(alpha, a, Nsim);  % Simulate positive jump size
        end

        for kk = 1:max_n
            V_n(:, kk) = simulateV(alpha, a, Nsim);  % Simulate negative jump size
        end

        beta_hat_p = beta_p * V_p;  % Adjusted positive beta
        beta_hat_n = beta_n * V_n;  % Adjusted negative beta

        j_hat_p = random('Gamma', (1 - alpha) * ones(size(beta_hat_p)), 1 ./ beta_hat_p);  % Positive Gamma increments
        j_hat_n = random('Gamma', (1 - alpha) * ones(size(beta_hat_n)), 1 ./ beta_hat_n);  % Negative Gamma increments

        count_matrix_p = repmat(1:max_p, Nsim, 1);  % Positive jump count matrix
        count_matrix_n = repmat(1:max_n, Nsim, 1);  % Negative jump count matrix

        rep_matrix_p = repmat(poisson_p(:, ii), 1, max_p);  % Repeated Poisson positive jump matrix
        rep_matrix_n = repmat(poisson_n(:, ii), 1, max_n);  % Repeated Poisson negative jump matrix

        find_p = rep_matrix_p >= count_matrix_p;  % Find valid positive jumps
        find_n = rep_matrix_n >= count_matrix_n;  % Find valid negative jumps

        % Computation of the theorical cumulants of the increment process
        % in order to compasate it 
        % Necessary step in order to join [4] Sabino and [1] Baviera &
        % Manzoni studies
        % Theoretical positive cumulants of the increment process
        cumulants_p = ctsCumulants(0, alpha, beta_p, c_p, dt, b, flag_OU_TS_FV);  
        % Theoretical negative cumulants of the increment process
        cumulants_n = ctsCumulants(0, alpha, beta_n, c_n, dt, b, flag_OU_TS_FV); 

        X2_p = sum(j_hat_p .* find_p, 2);  % Sum of positive increments
        X2_n = sum(j_hat_n .* find_n, 2);  % Sum of negative increments
        
        % Updating the process
        X(:, ii+1) = gamma_c * dt + a * X(:, ii) + X1_p + X2_p - cumulants_p(1) ...
                     - (X1_n + X2_n - cumulants_n(1));  
                 
        % Computing the drift term
        drift(ii) = -psi_Z_OU_TS(-1i, dt * ii, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c);  

    end

    %% Computing cumulants

    % Checking the last time step
    simCumulantsT = compCumulants(X(:, M+1)) * 1000; 
    theorCumulantsT = bctsCumulants(x0, alpha, beta_p, beta_n, c_p, c_n, gamma_c, T, b, flag_OU_TS_FV) * 1000;  

    % Checking the first time step
    simCumulants_dt = compCumulants(X(:, 2)) * 1000;  
    theorCumulants_dt = bctsCumulants(x0, alpha, beta_p, beta_n, c_p, c_n, gamma_c, dt, b, flag_OU_TS_FV) * 1000; 

    %% Computing the log forward in the risk neutral measure   
    logFwd = X + repmat([0, drift], Nsim, 1);  % Compute log forward prices

end % function sim_OU_TS_FinVar_ED
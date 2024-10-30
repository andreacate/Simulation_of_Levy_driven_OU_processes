function [X, theorCumulantsT, simCumulantsT, theorCumulants_dt, simCumulants_dt, logFwd] ...
                        = sim_OU_TS_FinAct_ED(x0, alpha, b, ...
                         beta_p, beta_n, c_p, c_n, gamma_c, T, M, Nsim)
% Simulation of OU-TS Finite Activity with Exact Decomposition following
% the algorithm 1 in Sabino [3]
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
% function theorCumulants(X0, alpha, beta, c, dt, b, flag)
% function bctsCumulants(X0, alpha, beta_p, beta_n, c_p, c_n, gamma_c, dt, b, flag)
% function compCumulants(vec)
% function psi_Z_OU_TS(u, dt, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c)

    %% Quantities of interest
    
    flag_OU_TS_FA = 1;  % Flag for OU-TS Finite Activity model
    X = zeros(Nsim, M+1);  % Initialize the matrix to store the simulations
    X(:,1) = x0*ones(Nsim,1);  % Set the initial condition for all simulations
    drift = zeros(1,M);  % Initialize the drift vector
    dt = T/M;  % Time step size

    % Lambda computation
    Gamma = gamma(-alpha);  
    a = exp(-b*dt);  
    lambda_p = c_p * Gamma * beta_p^(alpha);  
    lambda_n = c_n * Gamma * beta_n^(alpha);  

    %% Simulations
    
    % Poisson simulation
    poisson_p = random('Poisson', lambda_p*dt, Nsim, M);  % Positive Poisson process
    poisson_n = random('Poisson', lambda_n*dt, Nsim, M);  % Negative Poisson process
    
    % Process simulation
    for ii = 1:M

        % Compound Poisson simulation
        max_p = max(poisson_p(:,ii));  % Maximum number of positive jumps
        max_n = max(poisson_n(:,ii));  % Maximum number of negative jumps

        uniform_p = rand(Nsim, max_p);  % Uniform random variables for positive jumps
        uniform_n = rand(Nsim, max_n);  % Uniform random variables for negative jumps

        beta_hat_p = beta_p * exp(b*uniform_p*dt);  % Adjusted beta for positive jumps
        beta_hat_n = beta_n * exp(b*uniform_n*dt);  % Adjusted beta for negative jumps

        j_hat_p = random('Gamma', -alpha*ones(size(uniform_p)), 1./beta_hat_p);  % Gamma-distributed jumps for positive jumps
        j_hat_n = random('Gamma', -alpha*ones(size(uniform_n)), 1./beta_hat_n);  % Gamma-distributed jumps for negative jumps

        count_matrix_p = repmat([1:max_p], Nsim, 1);  % Matrix for counting positive jumps
        count_matrix_n = repmat([1:max_n], Nsim, 1);  % Matrix for counting negative jumps

        rep_matrix_p = repmat(poisson_p(:,ii), 1, max_p);  % Replicated Poisson process for positive jumps
        rep_matrix_n = repmat(poisson_n(:,ii), 1, max_n);  % Replicated Poisson process for negative jumps

        find_p = rep_matrix_p >= count_matrix_p;  % Logical matrix for valid positive jumps
        find_n = rep_matrix_n >= count_matrix_n;  % Logical matrix for valid negative jumps

        X2_p = sum(j_hat_p .* find_p, 2);  % Sum of positive jumps
        X2_n = sum(j_hat_n .* find_n, 2);  % Sum of negative jumps

        % Computation of the theorical cumulants of the increment process
        % in order to compasate it 
        % Necessary step in order to join [3] Sabino and [1] Baviera &
        % Manzoni studies
        % Theoretical positive cumulants of the increment process
        cumulants_p = ctsCumulants(0, alpha, beta_p, c_p, dt, b, flag_OU_TS_FA); 
        % Theoretical negative cumulants of the increment process
        cumulants_n = ctsCumulants(0, alpha, beta_n, c_n, dt, b, flag_OU_TS_FA);  
        
        % Updating the process
        X(:,ii+1) = gamma_c * dt + X(:,ii) * a + X2_p - cumulants_p(1) - (X2_n - cumulants_n(1));  
        % Compute drift term 
        drift(ii) = -psi_Z_OU_TS(-1i, dt*ii, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c);  
    end

    %% Computing cumulants

    % Checking the last time step
    simCumulantsT = compCumulants(X(:, M+1)) * 1000;  
    theorCumulantsT = bctsCumulants(x0, alpha, beta_p, beta_n, c_p, c_n, gamma_c, T, b, flag_OU_TS_FA) * 1000;  

    % Checking the first time step
    simCumulants_dt = compCumulants(X(:, 2)) * 1000; 
    theorCumulants_dt = bctsCumulants(x0, alpha, beta_p, beta_n, c_p, c_n, gamma_c, dt, b, flag_OU_TS_FA) * 1000;

    %% Computing the log forward in the risk neutral measure

    logFwd = X + repmat([0, drift], Nsim, 1);  

end % function sim_OU_TS_FinAct_ED
function [X, theorCumulantsT, simCumulantsT, theorCumulants_dt, simCumulants_dt, logFwd] = sim_Gaussian_OU_ED(x0, b, sigma, T, M, Nsim)
% Simulation for Gaussian OU process with exact method
%
% INPUT:
% x0:       initial condition
% b:        mean reverting parameter
% sigma:    volatility
% T:        time to maturity
% M:        number of steps
% Nsim:     number of simulations
% 
% OUTPUT:
% X:                   process simulated
% cumulants_process:   theoretical cumulants
% cumulants_sim:       simulated cumulants
% logFwd:              log of fwd prices
% 
% USES:
% function compCumulants(vec)

   
    %% Pre-allocate memory for simulations
    % Allocate memory for X (half for original, half for antithetic)
    X = zeros(Nsim/2, M+1);
    X(:, 1) = x0 * ones(Nsim/2, 1);
    
    % Allocate memory for antithetic process (X_AV)
    X_AV = zeros(Nsim/2, M+1);
    X_AV(:, 1) = x0 * ones(Nsim/2, 1);
    
    %% Simulation parameters
    % Time step for each simulation step
    dt = T/M;
    % Variance per diffusion step (considers mean reversion effect)
    var_dx = (exp(2*b*dt) -1 ) * (sigma^2*exp(-2*b*dt))/(2*b);
    
    %% Simulation loop
    % Loop through each time step
    for ii = 1:M
        % Generate standard normal random numbers
        Z = randn(Nsim/2,1);
        % Simulate X process' strong solution
        X(:,ii+1) = X(:,ii)*exp(-b*dt) + sqrt(var_dx) * Z;
        % Simulate antithetic X_AV process using same parameters but opposite noise (-Z)
        X_AV(:,ii+1) = X_AV(:,ii)*exp(-b*dt) - sqrt(var_dx) * Z;
    end
    
    %% Combine simulations and calculate drift
    % Combine simulated paths (X) and antithetic paths (X_AV)
    X = [X; X_AV];
    % Drift for each time step, correction to be in the RN measure
    drift = 0.5 *sigma^2/(2*b)*(1-exp(-2*b*[0:dt:T]));
    % Calculate log-forward prices 
    logFwd = X - drift;
    
    %% Cumulants
    % All scaled by 1000 for readability
    
    % Theoretical cumulants for the first step 
    theorCumulants_dt = [0; sigma^2/(2*b)*(1-exp(-2*b*dt)); 0; 0] * 1000;
    % Calculate simulated for the first step 
    simCumulants_dt = compCumulants(X(:,2)) * 1000;

    % Theoretical cumulants for the first step 
    theorCumulantsT = [0; sigma^2/(2*b)*(1-exp(-2*b*dt*M)); 0; 0] * 1000;
    % Calculate simulated for the last step 
    simCumulantsT = compCumulants(X(:, M+1)) * 1000;  
    

end % function sim_Gaussian_OU_ED

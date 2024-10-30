function plotSimulation (x0, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c, T, M, Nsim, model)
% Plot of the time evolution of the process simulations
%
% INPUT:
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
% model:    1 -> OU-CTS Finite Activity
%           2 -> CTS-OU Finite Variation
%           3 -> OU-CTS Finite Variation
% 
% USES
% function sim_OU_TS_FinAct_ED(x0, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c, T, M, Nsim)
% function sim_OU_TS_FinVar_ED(x0, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c, T, M, Nsim)
% function sim_TS_OU_FinVar_ED(x0, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c, T, M, Nsim)

    %% Generate simulations

    switch model
      case 1
        [X, ~, ~, ~, ~, ~] = sim_OU_TS_FinAct_ED(x0, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c, T, M, Nsim);
      case 2
        [X, ~, ~, ~, ~, ~] = sim_TS_OU_FinVar_ED(x0, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c, T, M, Nsim);
      case 3
        [X, ~, ~, ~, ~, ~] = sim_OU_TS_FinVar_ED(x0, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c, T, M, Nsim);
    end

    %% Generate plot

    figure() 
    grid on 

    % Plot all simulated paths (each column represents a single simulation)
    plot(linspace(0, T, M+1), X(:,:)) 

    switch model
      case 1
        title('OU-CTS Finite Activity Simulations')
      case 2
        title('TS-OU Finite Variation Simulations')
      case 3
        title('OU-CTS Finite Variation Simulations')
    end
    
end % function plotSimulation
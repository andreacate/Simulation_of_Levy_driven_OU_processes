function cumulants = bctsCumulants(X0, alpha, beta_p, beta_n, c_p, c_n, gamma_c, dt, b, flag)
% Bilateral TS Cumulants Computation for BCTS-OU and OU-BCTS 
% as described in [1] Baviera
%
% INPUT
% X0:      initial condition
% alpha:   stability parameter    
% beta_p:  positive beta
% beta_n:  negative beta 
% c_p:     positive c 
% c_n:     negative c 
% gamma_c: drift 
% dt:      time interval
% b:       mean reverting parameter
% flag:    1 -> OU-CTS Finite Activity
%          2 -> CTS-OU Finite Variation
%          3 -> OU-CTS Finite Variation

% NOTE: We added (-1)^k in every case

    %% Quantities of interest
    cumulants = zeros(4,1);  % Initialize the cumulants vector
    cumulants_L = zeros(4,1);  % Initialize the cumulants of L vector
    cumulants_X = zeros(4,1);  % Initialize the cumulants of X vector
    k = [1:4]';  % Vector of indices 1 to 4

    %% Cumulants computation
    
    switch flag  % Switch based on the flag value to choose the appropriate model

        case 1  % Case for OU-CTS Finite Activity

        cumulants_L(1) = gamma_c;  % Set the first cumulant of L to gamma_c

        % Compute the cumulants for k=2,3,4
        cumulants_L(2:end) = c_p*beta_p.^(alpha-k(2:end)).*gamma(k(2:end)-alpha) ...
                             + (-1).^k(2:end) .* c_n .* beta_n.^(alpha-k(2:end)).*gamma(k(2:end)-alpha);

        % Compute the overall cumulants
        cumulants = X0*exp(-b*dt).*(k == ones(size(k))) + cumulants_L./(b*k) .* (1-exp(-k*b*dt));

        case 2  % Case for CTS-OU Finite Variation

        cumulants_X(1) = gamma_c;  % Set the first cumulant of X to gamma_c

        % Compute the cumulants for k=2,3,4
        cumulants_X(2:end) = c_p*beta_p.^(alpha-k(2:end)).*gamma(k(2:end)-alpha) ...
                             + (-1).^k(2:end) .* c_n .* beta_n.^(alpha-k(2:end)).*gamma(k(2:end)-alpha);

        % Compute the overall cumulants
        cumulants = X0*exp(-b*dt).*(k == ones(size(k))) + cumulants_X .* (1-exp(-k*b*dt));

        % Alternative version commented out
        % cumulants = X0*exp(-b*dt).*(k == ones(size(k))) + c_p * beta_p.^(alpha-k) .* gamma(k-alpha) .* (1-exp(-k*b*dt))...
        %                      + (-1).^k .* c_n .* beta_n.^(alpha-k) .* gamma(k-alpha) .* (1-exp(-k*b*dt));
        % cumulants(1)=gamma_c;

        case 3  % Case for OU-CTS Finite Variation

        cumulants_L(1) = gamma_c;  % Set the first cumulant of L to gamma_c
        % Compute the cumulants for k=2,3,4
        
        cumulants_L(2:end) = c_p*beta_p.^(alpha-k(2:end)).*gamma(k(2:end)-alpha) ...
                             + (-1).^k(2:end) .* c_n .* beta_n.^(alpha-k(2:end)).*gamma(k(2:end)-alpha);

        % Compute the overall cumulants
        cumulants = X0*exp(-b*dt).*(k == ones(size(k))) + cumulants_L./(b*k) .* (1-exp(-k*b*dt));
        
    end  
 
end % End of function bctsCumulants
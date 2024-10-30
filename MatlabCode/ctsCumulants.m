function cumulants = ctsCumulants(X0, alpha, beta, c, dt, b, flag)
% Theoretical cumulants computation for CTS-OU and OU-CTS 
% as described in [3] Sabino, [4] Sabino.
%
% INPUT
% X0:      initial condition
% alpha:   stability parameter    
% beta:    beta 
% c:       model parameter
% dt:      time interval
% b:       mean reverting parameter
% flag:    1 -> Finite Activity
%          2 -> CTS-OU Finite Variation
%          3 -> OU-CTS Finite Variation

    %% Quantities of interest
    cumulants = zeros(4,1);
    k = [1:4]';
    
    %% Cumulants computation
    
    switch flag

        case 1
        cumulants_L = c * beta.^(alpha-k) .* gamma(k-alpha);
        cumulants = X0*exp(-b*dt).*(k == ones(size(k))) + cumulants_L./(b*k) .* (1-exp(-k*b*dt));

        case 2
        cumulants_X = c * beta.^(alpha-k) .* gamma(k-alpha);
        cumulants = X0*exp(-b*dt).*(k == ones(size(k))) + cumulants_X .* (1-exp(-k*b*dt));

        case 3
        % We added b to the exponential
        cumulants_L = c * beta.^(alpha-k) .* gamma(k-alpha);
        cumulants = X0*exp(-b*dt).*(k == ones(size(k))) + cumulants_L./(b*k) .* (1-exp(-k*b*dt));       

    end
    
end % function theorCumulants
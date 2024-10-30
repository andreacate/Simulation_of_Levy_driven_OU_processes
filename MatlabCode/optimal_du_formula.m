function du = optimal_du_formula(alpha, a, b, c_p, c_n, dt, M, scale, flag)
% Optimal du by using the formula in Baviera [1] based on making of comparable 
% order of the truncation and discretization error
% 
% INPUT:
% alpha:    stability parameter    
% a:        shift
% b:        mean reverting parameter
% beta_p:   positive beta
% beta_n:   negative beta
% c_p:      positive c
% c_n:      negative c
% dt:       time step
% M:        number of steps
% scale:    scaling factor
% flag:     1 -> OU-CTS Finite Activity
%           2 -> CTS-OU Finite Variation
%           3 -> OU-CTS Finite Variation
%

    %% Quantities of interest
    
    N = 2^M;
    
    
    %% Computation of du
    
    if flag == 1
        % This block tries to calculate du for the OU-TS Finite Activity model (commented out).
        % The commented code uses an optimization routine to find du.
        % An alternative fixed value (e.g., du = 0.4) is provided but might not be optimal.
        %             % Lambda computation
        %             Gamma = gamma(-alpha);
        %             a = exp(-b*dt);
        %             lambda_p = c_p * Gamma * beta_p^(alpha);
        %             lambda_n = c_n * Gamma * beta_n^(alpha);
        %
        %             B=(c_p/(exp(lambda_p*dt)-1)+c_n/(exp(lambda_n*dt)-1))...
        %                 *Gamma*(1-exp(-alpha*b*dt))/(alpha*b);
        %
        %             f = @(h) (N*h)^alpha-exp(2*pi*abs(a)/h);
        %             options = optimoptions('lsqnonlin', ...
        %                                'TolFun', 1e-8);
        %
        %             du = lsqnonlin(f, [0.2], [0.001], [0.5]);
    
    elseif flag == 2
        % optimal du calculation for the TS-OU Finite Variation model.
        l = -(c_p+c_n) * gamma(-alpha) * cos(alpha*pi/2) * (1-exp(-alpha*b*dt)) * scale^alpha;
        du = ((2*pi*abs(a))/(l*N^alpha))^(1/(alpha+1));
    elseif flag == 3
        % optimal du calculation for the OU-TS Finite Variation model.
        l = -(c_p+c_n) * gamma(-alpha) * cos(alpha*pi/2) * (1-exp(-alpha*b*dt))/(alpha*b) * scale^alpha;
        du = ((2*pi*abs(a))/(l*N^alpha))^(1/(alpha+1));
    end

end % function optimal_du_formula
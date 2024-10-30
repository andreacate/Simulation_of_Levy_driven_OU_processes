function psi = psi_V(u, dt, alpha, b, beta_p, beta_n, c_p, c_n)
% Characteristic exponent of process V for FA processes ([1] Baviera)
% 
% INPUT: 
% u:         vector in which we compute psi V
% dt:        time step
% alpha:     stability parameter    
% b:         mean reverting parameter
% beta_p:    positive beta
% beta_n:    negative beta
% c_p:       positive c
% c_n:       negative c

    % Lambda computation - Calculate the total and individual jump intensities (lambda)
    Gamma = gamma(-alpha);  
    lambda_p = c_p * Gamma * beta_p^(alpha);  
    lambda_n = c_n * Gamma * beta_n^(alpha);  
    lambda_tot = lambda_p + lambda_n;  
    
    % Characteristic function - Calculate the characteristic functions for positive and negative jumps
    phi_J_p = 1/(b*dt) * integral(@(z) (z-1i*u).^alpha./z.^(alpha+1), beta_p, beta_p*exp(b*dt),"ArrayValued",true); 
    phi_J_n = 1/(b*dt) * integral(@(z) (z+1i*u).^alpha./z.^(alpha+1), beta_n, beta_n*exp(b*dt),"ArrayValued",true); 
    
    % Psi calculation - Combine jump intensities and characteristic functions to compute psi
    psi = log( (exp((lambda_tot*dt) .* (1/(lambda_tot) .* (lambda_p.*phi_J_p  + lambda_n.*phi_J_n))) - 1) ./ ...
             (exp(lambda_tot*dt)-1) );
    
end % function psi_V
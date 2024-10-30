function psi = psi_Z_OU_TS(u, dt, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c)
% Characteristic exponent of process Z in OU-TS processes 
% ([1] Baviera & Manzoni)
% 
% INPUT:
% u:         vector in which we compute psi Z
% alpha:     stability parameter  
% b:         mean reverting parameter
% beta_p:    positive beta
% beta_n:    negative beta
% c_p:       positive c
% c_n:       negative c
% gamma_c:   drift

    %% Characteristic exponent
    
    int_p = integral(@(z) (z-1i*u).^alpha./z.^(alpha+1), beta_p, beta_p*exp(b*dt), "ArrayValued", true);
    int_n = integral(@(z) (z+1i*u).^alpha./z.^(alpha+1), beta_n, beta_n*exp(b*dt), "ArrayValued", true);

    psi = 1i*u.*(1-exp(-b*dt))/b * gamma_c ...
          + c_p *beta_p^alpha*gamma(-alpha)/b .* ...
          (int_p - b*dt + alpha/beta_p *1i*u* (1-exp(-b*dt))) ...
          + c_n * beta_n^alpha*gamma(-alpha)/b .* ...
          (int_n - b*dt - alpha/beta_n *1i*u* (1-exp(-b*dt)));

end % function psi_Z_OU_TS
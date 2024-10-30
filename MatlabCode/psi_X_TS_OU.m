function psi = psi_X_TS_OU(u, alpha, beta_p, beta_n, c_p, c_n, gamma_c)
% Characteristic exponent of process X in TS-OU processes 
% ([1] Baviera & Manzoni)
% 
% INPUT:
% u:         vector in which we compute psi V
% alpha:     stability parameter    
% beta_p:    positive beta
% beta_n:    negative beta
% c_p:       positive c
% c_n:       negative c
% gamma_c:   drift

    %% Characteristic exponent
    
    psi = 1i*u*gamma_c + c_p * gamma(-alpha) * beta_p^alpha * ((1-1i*u/beta_p).^alpha -1+1i*u*alpha/beta_p) ...
          + c_n * gamma(-alpha) * beta_n^alpha * ((1+1i*u/beta_n).^alpha -1 -1i*u*alpha/beta_n);
         

end % function psi_X_TS_OU
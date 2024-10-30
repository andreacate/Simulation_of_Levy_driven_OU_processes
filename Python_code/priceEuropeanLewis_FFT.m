function prices = priceEuropeanLewis_FFT(S0, b, alpha, beta_p, beta_n, c_p, c_n, gamma_c, T, moneyness, r, flag, M_fft)
% Computation of European Call Price via Lewis formula
%
% INPUT:
% S0:        initial underlying condition  
% b:         mean reverting parameter
% alpha:     stability parameter  
% beta_p:    positive beta
% beta_n:    negative beta
% c_p:       positive c
% c_n:       negative c
% gamma_c:   drift
% T:         time to maturity
% moneyness: vector of moneyness
% r:         rate
% flag:      1 -> OU-CTS Finite Activity
%            2 -> CTS-OU Finite Variation
%            3 -> OU-CTS Finite Variation
% M_fft:     parameter of FFT
%
% USES:
% function psi_Z_OU_TS(u, dt, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c)
% function psi_X_TS_OU(u, alpha, beta_p, beta_n, c_p, c_n, gamma_c)
% function integralViaFFT(f, M, x1)

    % Integrand computation - Define integrand based on chosen model (flag)
    
    switch flag
        case 1
            % For OU-CTS Finite Activity model
            drift = -psi_Z_OU_TS(-1i,T,alpha,b,beta_p,beta_n,c_p,c_n,gamma_c);
            integrand = @(u) exp(psi_Z_OU_TS(-u-1i/2,T,alpha,b,beta_p,beta_n,c_p,c_n,gamma_c))./(u.^2+1/4).*exp(1i*(-u-1i/2).*drift);
        case 2
            % For CTS-OU Finite Variation model
            drift = -(psi_X_TS_OU(-1i,alpha,beta_p,beta_n,c_p,c_n,gamma_c)-psi_X_TS_OU(-1i*exp(-b*T),alpha,beta_p,beta_n,c_p,c_n,gamma_c));
            integrand = @(u) exp(psi_X_TS_OU(-u-1i/2,alpha,beta_p,beta_n,c_p,c_n,gamma_c) - ...
                             psi_X_TS_OU((-u-1i/2)*exp(-b*T),alpha,beta_p,beta_n,c_p,c_n,gamma_c))./...
                             (u.^2+1/4).*exp(1i*(-u-1i/2).*drift);
        case 3
            % For OU-CTS Finite Variation model (same as case 1)
            drift = -psi_Z_OU_TS(-1i,T,alpha,b,beta_p,beta_n,c_p,c_n,gamma_c);
            integrand = @(u) exp(psi_Z_OU_TS(-u-1i/2,T,alpha,b,beta_p,beta_n,c_p,c_n,gamma_c))./(u.^2+1/4).*exp(1i*(-u-1i/2).*drift);
    end
    
    % Integral via FFT - Perform numerical integration using FFT  
    [integral_fft, x_grid] = integralViaFFT_x1(integrand, M_fft, moneyness(1)*1000);
    integral_Lewis = interp1(x_grid, integral_fft, moneyness, 'spline');
    
    % European price - Calculate European call option price using Lewis formula
    fwd = S0 * exp(r*T);
    prices = exp(-r*T) * fwd * (1 - exp(-moneyness/2) .* real(integral_Lewis)/(2*pi));
                          
end % function priceEuropeanLewis_FFT
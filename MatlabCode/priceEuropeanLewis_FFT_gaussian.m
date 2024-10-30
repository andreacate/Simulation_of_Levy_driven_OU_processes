function price = priceEuropeanLewis_FFT_gaussian(S0, b, sigma, T, moneyness, r,  M_fft)
% Computation of European Call Price via Lewis formula for Gaussian OU
% process
%
% INPUT:
% S0:        initial underlying condition  
% b:         mean reverting parameter
% sigma:     volatility
% T:         time to maturity
% moneyness: vector of moneyness
% r:         rate
% M_fft:     parameter of FFT
%
% USES:
% function integralViaFFT(f, M, x1)

    %% Quantities of interest
    
    % Variance per diffusion step (considers mean reversion effect)
    var_dx = (exp(2*b*T) -1 )* (sigma^2*exp(-2*b*T))/(2*b);
    
    % Characteristic exponent function (psi_Z) for Gaussian OU
    psi_Z = @(u)  - (u).^2*var_dx/2;
    
    % Drift calculation using psi_Z evaluated at -1i
    drift = -psi_Z(-1i);

    % Alternative formula for the drift, equivalent under Gaussian OU (commented out)
    % drift = 0.5 *sigma^2 * (1-exp(-2*b*T)) / (2*b);
    
    %% Integrand definition
    
    % Integrand function used in the Lewis formula
    integrand=@(u) exp(psi_Z(-u-1i/2))./(u.^2+1/4).*exp(1i*(-u-1i/2).*drift);
   
    %% Integral via FFT
    % Perform numerical integration using the integralViaFFT function
    [integral_fft, x_grid] = integralViaFFT_x1(integrand, M_fft, moneyness(1)*1500);
    % Interpolation for the Lewis integral using spline interpolation
    integral_Lewis = interp1(x_grid, integral_fft, moneyness, 'spline');
    
    %% European price calculation
    % Forward price considering risk-free interest rate
    fwd = S0 * exp(r*T);
    % European call price using Lewis formula (real part of integral_Lewis)
    price = exp(-r*T) * fwd * (1 - exp(-moneyness/2) .* real(integral_Lewis)/(2*pi));
    
end % function priceEuropeanLewis_FFT_gaussian
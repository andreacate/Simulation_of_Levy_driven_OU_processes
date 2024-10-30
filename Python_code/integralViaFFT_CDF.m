function [I, x] = integralViaFFT_CDF(M_fft, du, dt, alpha, a, b, beta_p,...
                                     beta_n, c_p, c_n, gamma_c, scale, flag)
% Function which computes the integral in the formula for the CDF in
% Baviera [1] using the FFT procedure 
% 
% INPUT:
% M_fft:        power of 2 to compute the number of intervals  
% du:           step of the integral grid
% dt:           time step
% alpha:        stability parameter  
% a:            shift
% beta_p:       positive beta
% beta_n:       negative beta
% c_p:          positive c
% c_n:          negative c
% gamma_c:      drift
% scale:        scaling factor
% flag:         1 -> OU-CTS Finite Activity
%               2 -> CTS-OU Finite Variation
%               3 -> OU-CTS Finite Variation
%
% OUTPUT
% I:            integral with FFT
% x:            x grid
%
% USES
% function psi_V(u, dt, alpha, b, beta_p, beta_n, c_p, c_n)
% function psi_X_TS_OU(u, alpha, beta_p, beta_n, c_p, c_n, gamma_c)
% function psi_Z_OU_TS(u, dt, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c)


%% FFT parameters:
% Number of data points used in the FFT calculation (power of 2)
N = 2^M_fft;    
% Spacing between the input of the CDF
dx = 2*pi/(du*N);   
% Starting point for the input of the CDF
x1 = -(N-1)*dx/2;  
% Array of data points ranging from x1 to its negative counterpart
x = [x1:dx:-x1];    
% Starting point for the frequency domain
u1 = -du*(N-1)/2; 
% Array of integral domain
u = [u1:du:-u1];   


%% Integrand computation:
% This step defines the integrand based on the chosen model (specified by flag).
% It utilizes helper functions (psi_V, psi_X_TS_OU, psi_Z_OU_TS) to calculate the characteristic functions.

switch flag
  case 1
    % OU-TS Finite Activity model
    integrand = exp(psi_V(scale*(u+1i*a), dt, alpha, b, beta_p, beta_n, c_p, c_n)) ./ (1i*(u+1i*a));
  case 2
    % TS-OU Finite Variation model
    psi_Z = (psi_X_TS_OU(scale*(u+1i*a), alpha, beta_p, beta_n, c_p, c_n, gamma_c) - ...
             psi_X_TS_OU(scale*(u+1i*a)*exp(-b*dt), alpha, beta_p, beta_n, c_p, c_n, gamma_c));
    integrand = exp(psi_Z) ./ (1i*(u+1i*a));
  case 3
    % OU-TS Finite Variation model
    integrand = exp(psi_Z_OU_TS(scale*(u+1i*a), dt, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c)) ./ (1i*(u+1i*a));
end

%% Integral computation:
% This step performs the FFT calculation and multiplies the result with a prefactor 
% to obtain the integral approximation. The real part is extracted as the integral value 

% Compute the Fast Fourier Transform:
FFT = fft(integrand.*exp(-1i.*x1.*du*[0:N-1]));

% Compute the integral by multiplying for the prefactor:
I = real(du * exp( -1i .* u1 .* x) .* FFT);

end % function integralViaFFT_CDF
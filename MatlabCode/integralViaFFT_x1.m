function [I, x] = integralViaFFT_x1(f, M, x1)
% Function which computes the integral in the Lewis formula for the price
% of a call option via the Fast Fourier Transform given the initial point
% x1
% 
% INPUTS:
% f:            function to be Fourier transformed         
% M:            power of 2 to compute the number of intervals          
% x1:           first value in the function domain
% 
% OUTPUT:
% I:            value of the integral

    %% FFT parameters:
    N = 2^M;    % points in the fft discretization
    dx = -2*x1/(N-1);   % interval in function domain
    x = [x1:dx:-x1];    % function domain discretization
    du = 2*pi/(dx*N);   % step for the integration
    u1 = -du*(N-1)/2;    % first point in integration domain
    u = [u1:du:-u1];     % integration domain discretization

    %% Compute the the Fast Fourier Transform:
    FFT = fft(f(u).*exp(-1i.*x1.*du*[0:N-1]));  % fft computation after phase correction

    %% Compute the integral by multiplying for the prefactor:
    I = real(du * exp( -1i .* u1 .* x) .* FFT); % prefactor multiplication for Fourier transform computation

end % function integralViaFFT





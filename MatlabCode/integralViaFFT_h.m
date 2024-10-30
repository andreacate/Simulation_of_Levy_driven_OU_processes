function [I, x] = integralViaFFT_h(f,M,du)
% Function which computes the integral in the Lewis formula for the price
% of a call option via the Fast Fourier Transform given the integral 
% discretization step du
% 
% INPUTS:
% f:            function to be Fourier transformed         
% M:            power of 2 to compute the number of intervals          
% dx:           step of the moneyness grid
% OUTPUT:
% I:            value of the integral

    %% FFT parameters:
    N = 2^M;
    dx = 2*pi/(du*N);
    x1 = -(N-1)*dx/2;
    x = [x1:dx:-x1];
    u1 = -du*(N-1)/2;
    u = [u1:du:-u1];

    %% Integral computation
    
    % Compute the the Fast Fourier Transform:
    FFT = fft(f(u).* exp(-1i.* x1.* du* [0:N-1] ));
    % Compute the integral by multiplying for the prefactor:
    I = real(du * exp( -1i .* u1 .* x) .* FFT);

end % function integralViaFFT_h
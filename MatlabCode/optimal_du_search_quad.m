function [du_optimal] = optimal_du_search_quad(M_fft, dt, alpha, a, b, beta_p, beta_n, c_p, c_n, gamma_c, x_values_quad, integral_quad, flag)
% Search of the optimal du to use for the FFT by using a comparison between
% the integral computed via quadrature and via FFT
%
%
% INPUT:
% M_fft:        power of 2 to compute the number of intervals  
% dt:           time step
% alpha:        stability parameter  
% a:            shift
% beta_p:       positive beta
% beta_n:       negative beta
% c_p:          positive c
% c_n:          negative c
% gamma_c:      drift
% x_values_quad grid of integration for the quadrature
% integral_quad quadrature integral values 
% flag:         1 -> OU-CTS Finite Activity
%               2 -> CTS-OU Finite Variation
%               3 -> OU-CTS Finite Variation
%

    % du possible values 
    du_values = [0.001:0.02:0.8];

    % Initializing the vector of error
    err = zeros(size(du_values));

    % Computation of the different function to be integrate via FFT depending
    % on the model 
    if flag==1
        f_FFT = @(u) exp( psi_V(u+1i*a, dt, alpha, b, beta_p, beta_n, c_p, c_n)) ./ (1i*(u+1i*a));
    elseif flag==2
        f_FFT = @(u) exp( psi_X_TS_OU((u+1i*a), alpha, beta_p, beta_n, c_p, c_n, gamma_c) - ...
                    psi_X_TS_OU((u+1i*a)*exp(-b*dt), alpha, beta_p, beta_n, c_p, c_n, gamma_c)) ./ (1i*(u+1i*a));
    elseif flag==3
        f_FFT = @(u) exp( psi_Z_OU_TS((u+1i*a), dt, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c)) ./ (1i*(u+1i*a));
    end

    % Computation of the FFT integral for each possible values of du 
    for ii = 1:length(du_values)
        [integral_FFT_u, x_grid] = integralViaFFT_h(f_FFT ,M_fft,du_values(ii));
        % Filling the vector of the error between the quadrature integral and
        % the FFT integral (we need to use spline interpolation since the grid
        % are different)
        err(ii) = norm(interp1(x_grid, integral_FFT_u, x_values_quad, 'spline')-integral_quad);
    end


    figure()
    plot (du_values, err)
    legend('Error between FFT CDF and Quadrature CDF')
    title('Comparison between FFT CDF and Quadrature CDF')

    % Finding the index of the smaller error
    [~,index]=min(err);

    % Assing the optimal du between the one in the possible values
    du_optimal=du_values(index);

end % function optimal_du_search_quad

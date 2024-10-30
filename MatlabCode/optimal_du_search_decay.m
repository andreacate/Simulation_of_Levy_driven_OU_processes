function [du] = optimal_du_search_decay(a, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c, dt, M_fft, flag)
% Search of the optimal du to use for the FFT by using a comparison between
% the decay of the characteristic function and the theorical exponential decay using
% lsqnonlin function, as described in [1] Baviera.
%
% INPUT
% a:        vertical shift 
% alpha:    stability parameter  
% b:        mean reverting parameter
% beta_p:   positive beta
% beta_n:   negative beta 
% c_p:      positive c 
% c_n:      negative c 
% gamma_c:  drift 
% dt:       time interval
% M_fft:    parameter of FTT which controls the number of interval of the
%           grid
% flag:     1 -> OU-CTS Finite Activity
%           2 -> CTS-OU Finite Variation
%           3 -> OU-CTS Finite Variation

    % Generate logarithmically spaced values for u
    u_values = 10.^linspace(1, 10, 20);

    if flag == 1
        % OU-TS finite activity model
        % param:    
        %           1 - w
        %           2 - B
        % Define the function to minimize
        f = @(param) abs(abs(exp(psi_Z_OU_TS(u_values + 1i * a, dt, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c))) ...
            - param(2) * abs(u_values).^(-param(1)));
        % Optimize parameters using least squares non-linear fitting
        optimal_params = lsqnonlin(f, [0.2, 0.2], [0, 0], [2, 10]);

        % Plot the characteristic function decay and the theoretical decay
        figure()
        plot(u_values, abs(exp(psi_Z_OU_TS(u_values + 1i * a, dt, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c))))
        hold on
        plot(u_values, optimal_params(2) * u_values.^(-optimal_params(1)));
        legend('CF decay', 'Theoretical decay')

        % Set a fixed value for du
        du = 0.012;

    elseif flag == 2
        % TS-OU finite variation model
        % Compute the characteristic exponent difference
        psi_Z_TS_OU = psi_X_TS_OU(u_values + 1i * a, alpha, beta_p, beta_n, c_p, c_n, gamma_c) ...
            - psi_X_TS_OU((u_values + 1i * a) * exp(-b * dt), alpha, beta_p, beta_n, c_p, c_n, gamma_c);
        % param:    
        %           1 - l
        %           2 - w
        %           3 - B
        % Define the function to minimize
        f = @(param) abs(abs(exp(psi_Z_TS_OU)) ...
            - param(3) * exp(-param(1) * abs(u_values).^param(2)));
        % Optimize parameters using least squares non-linear fitting
        optimal_params = lsqnonlin(f, [0.2, 0.2, 0.2], [0, 0, 0], [10, 3, 100]);

        % Plot the characteristic function decay and the theoretical decay
        figure()
        semilogy(u_values, abs(exp(psi_Z_OU_TS(u_values + 1i * a, dt, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c))))
        hold on
        semilogy(u_values, optimal_params(3) * exp(-optimal_params(1) * u_values.^optimal_params(2)));
        legend('CF decay', 'Theoretical decay')

        % Calculate du based on the optimized parameters
        du = ((2 * pi * abs(a)) / (optimal_params(1) * (2^M_fft)^optimal_params(2)))^(1 / (optimal_params(2) + 1));

    elseif flag == 3
        % OU-TS finite variation model
        % param:    
        %           1 - l
        %           2 - w
        %           3 - B
        % Define the function to minimize
        f = @(param) abs(abs(exp(psi_Z_OU_TS(u_values + 1i * a, dt, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c))) ...
            - param(3) * exp(-param(1) * abs(u_values).^param(2)));
        % Optimize parameters using least squares non-linear fitting
        optimal_params = lsqnonlin(f, [0.2, 0.2, 0.2], [0, 0, 0], [10, 3, 100]);

        % Plot the characteristic function decay and the theoretical decay
        figure()
        semilogy(u_values, abs(exp(psi_Z_OU_TS(u_values + 1i * a, dt, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c))))
        hold on
        semilogy(u_values, optimal_params(3) * exp(-optimal_params(1) * u_values.^optimal_params(2)));
        legend('CF decay', 'Theoretical decay')

        % Calculate du based on the optimized parameters
        du = ((2 * pi * abs(a)) / (optimal_params(1) * (2^M_fft)^optimal_params(2)))^(1 / (optimal_params(2) + 1));

    end

end % function optimal_du_search_decay
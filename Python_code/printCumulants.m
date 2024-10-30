function printCumulants (theorCumulants_T, simCumulants_T, theorCumulants_dt, simCumulants_dt, T, M)
% Print of theoretical and empirical cumulants
%
% INPUT:
% theorCumulants_T:     theoretical cumulants last step
% simCumulants_T:       empirical cumulants last step
% theorCumulants_dt:    theoretical cumulants first step
% simCumulants_dt:      empirical cumulants first step
% T:                    time step
% M:                    number of step

    %% Print
    % Cumulants last step
    fprintf("\n Statistical cumulants for the first step with delta t %.2f : \n",T/M);
    fprintf("Statistical cumulants order 1 = %.5f\n", simCumulants_dt(1));
    fprintf("Statistical cumulants order 2 = %.5f\n", simCumulants_dt(2));
    fprintf("Statistical cumulants order 3 = %.5f\n", simCumulants_dt(3));
    fprintf("Statistical cumulants order 4 = %.5f\n", simCumulants_dt(4));
    fprintf("\n Theorical cumulants for the first step with delta t %.2f : \n",T/M);
    fprintf("Theorical cumulants order 1 = %.5f\n", theorCumulants_dt(1));
    fprintf("Theorical cumulants order 2 = %.5f\n", theorCumulants_dt(2));
    fprintf("Theorical cumulants order 3 = %.5f\n", theorCumulants_dt(3));
    fprintf("Theorical cumulants order 4 = %.5f\n", theorCumulants_dt(4));
    
    % Cumulants last step
    fprintf("\n Statistical cumulants for the last step with T %.2f : \n",T);
    fprintf("Statistical cumulants order 1 = %.5f\n", simCumulants_T(1));
    fprintf("Statistical cumulants order 2 = %.5f\n", simCumulants_T(2));
    fprintf("Statistical cumulants order 3 = %.5f\n", simCumulants_T(3));
    fprintf("Statistical cumulants order 4 = %.5f\n", simCumulants_T(4));
    fprintf("\n Theorical cumulants for the last step with T %.2f : \n",T);
    fprintf("Theorical cumulants order 1 = %.5f\n", theorCumulants_T(1));
    fprintf("Theorical cumulants order 2 = %.5f\n", theorCumulants_T(2));
    fprintf("Theorical cumulants order 3 = %.5f\n", theorCumulants_T(3));
    fprintf("Theorical cumulants order 4 = %.5f\n", theorCumulants_T(4));
    
end % function printCumulants
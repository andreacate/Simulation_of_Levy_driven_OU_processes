function cumulants = compCumulants(vec)
% Empirical cumulants computation
%
% INPUT
% vec:    column vector of simulations

  %% Quantities of interest
    cumulants = zeros(4,1);  % Initialize the cumulants vector to store four cumulants

    %% Moments computation

    % The first cumulant is the mean of the vector
    cumulants(1) = mean(vec); 
    % The second cumulant is the second moment (variance)
    cumulants(2) = moment(vec, 2);  
    % The third cumulant is the third moment
    cumulants(3) = moment(vec, 3);  
    % The fourth cumulant is the fourth moment minus three times the square of the second moment
    cumulants(4) = moment(vec, 4) - 3 * moment(vec, 2)^2;  

end % function compCumulants

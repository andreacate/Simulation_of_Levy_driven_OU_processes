function [x_CDF, values_CDF] = selectCDF(x, discreteCDF, toll)
% Selection of CDF values that respect the theoretical conditions
% 
% INPUT
% x:            x grid
% discreteCDF:  values of discrete CDF 
% toll:         tollerance
% 
% OUTPUT
% x_CDF:        selected values of x
% values_CDF:   selected values of CDF
% Check if the difference between consecutive CDF values is greater than the tolerance 
% and if the values are within the range of 0 to 1 (inclusive)

    check = (discreteCDF(2:end) > (discreteCDF(1:end-1)) - toll) .* ...
            (discreteCDF(2:end) >= 0 ) .* ...
            (discreteCDF(2:end) <= 1);

    % Managing critical cases (when the first or the last point are in the longest
    % subsequence)
    check = [0 check 0];

    %% Finding invalid indices

    % Find the indices where the check vector is 0, indicating a violation of the conditions
    indexes = find(check == 0);

    %% Finding the longest valid interval

    % Calculate the difference between consecutive indices to find the longest valid interval
    [~, indexI] = max(indexes(2:end) - indexes(1:end-1));

    % Extract the valid x and CDF values within the longest valid interval
    x_CDF = x(indexes(indexI)+1 : indexes(indexI+1)-1);
    values_CDF = discreteCDF(indexes(indexI)+1 : indexes(indexI+1)-1);

end % function selectCDF
function [sampleValues] = samplingFromDiscreteCDF(discrete_x, discrete_CDF, N_sample)
% Sampling values from discrete CDF
% 
% INPUT
% discrete_x:           grid of X
% discrete_CDF:         values of discrete CDF
% N_sample:             number of samples

    %% Quantities of interest

    % Estimate parameters for the negative exponential tail of the CDF (leftmost part)
    a_n = discrete_CDF(1);  
    b_n = log(discrete_CDF(2)/discrete_CDF(1))/(discrete_x(2)-discrete_x(1));  

    % Estimate parameters for the positive exponential tail of the CDF (rightmost part)
    a_p = discrete_CDF(end);  
    b_p = log((1 - discrete_CDF(end-1))/(1-discrete_CDF(end)))/(discrete_x(end-1)-discrete_x(end)); 

    % Generate N_sample uniform random variables between 0 and 1
    u = rand(N_sample, 1);

    %% Inverting CDF and sampling

    % Remove duplicates from discrete_CDF (in case they exist) and keep the corresponding indices in ia
    [discrete_CDF_unique, ia] = unique(discrete_CDF);

    % Define an anonymous function for inverting the CDF using linear interpolation
    % This function considers three regions:
    %   - Linear interpolation within the main body of the CDF (between discrete_CDF(1) and discrete_CDF(end))
    %   - Negative exponential tail for z <= discrete_CDF(1)
    %   - Positive exponential tail for z >= discrete_CDF(end)
    inv_CDF = @(z) (interp1(discrete_CDF_unique, discrete_x(ia), z,'linear','extrap')) .* (z<=discrete_CDF(end)) .* (z>=discrete_CDF(1)) ...
                      + (log(z/a_n)/b_n+discrete_x(1)) .* (z<=discrete_CDF(1)) ...
                      + (log(-(z-1)/a_p)/b_p+discrete_x(end)) .* (z>=discrete_CDF(end));

    % You can uncomment the following lines to visualize the inverted CDF
    % plot(linspace(0,1,1000),inv_CDF(linspace(0,1,1000)))
    % figure()
    % plot(discrete_CDF,inv_CDF(discrete_CDF),'*')

    %% Sample values with antithetic variables

    % Generate samples using the inverted CDF function
    sampleValues_1 = inv_CDF(u);

    % Generate antithetic samples by reflecting the uniform random variables across 0.5
    sampleValues_2 = inv_CDF(1-u);

    % Combine the original samples and their antithetic counterparts
    sampleValues = [sampleValues_1; sampleValues_2];
    
end % function samplingFromDiscreteCDF
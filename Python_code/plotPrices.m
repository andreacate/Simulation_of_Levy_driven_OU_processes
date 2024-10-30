function plotPrices(moneyness, pricesLewis, prices_ED, conf_int_ED, prices_FM, conf_int_FM)
% Plots of prices with different method
% 
% INPUT:
% moneyness:       moneyness
% pricesLewis:     prices computed with Lewis formula
% prices_ED:       prices with Exact Decomposition
% conf_int_ED:     confidence intervals with Exact Decomposition
% prices_FM:       prices with Fast General Monte Carlo
% conf_int_FM:     confidence intervals with Fast General Monte Carlo

    figure()  

    % Plot the mid derivative prices from the Exact Decomposition method
    plot(moneyness, prices_ED)
    hold on 
    grid on  

    % Plot the confidence interval bounds for the Exact Decomposition method
    plot(moneyness, conf_int_ED(1,:), '--') 
    plot(moneyness, conf_int_ED(2,:), '--')  

    % Plot the mid derivative prices from the Fast General Monte Carlo method
    plot(moneyness, prices_FM)

    % Plot the confidence interval bounds for the Fast General Monte Carlo method
    plot(moneyness, conf_int_FM(1,:), '--')  
    plot(moneyness, conf_int_FM(2,:), '--')  

    % Plot the prices calculated using the Lewis formula 
    plot(moneyness, pricesLewis, '*-') 

    legend('Mid Derivative Prices ED', 'Lower Bound ED', 'Upper Bound ED',...
           'Mid Derivative Prices FM', 'Lower Bound FM', 'Upper Bound FM', 'Prices with Lewis', 'FontSize', 15)

end % function plotPrices

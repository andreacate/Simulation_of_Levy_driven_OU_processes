function plotVolatilities (pricesLewis, prices_ED, prices_FM, conf_int_ED, conf_int_FM, moneyness, r, T, S0)
% Plot of volatility surfaces
%
% INPUT:
% pricesLewis:     prices computed with Lewis formula
% prices_ED:       prices with Exact Decomposition
% prices_FM:       prices with Fast General Monte Carlo
% conf_int_ED:     confidence intervals with Exact Decomposition
% conf_int_FM:     confidence intervals with Fast General Monte Carlo
% moneyness:       moneyness
% r:               rate
% T:               time to maturity
% S0:              initial underlying condition
% 
% USES:
% plotVolSmile(prices,moneyness,rate,TTM,S0)

    %% Generate plots

    figure() 
    hold on
    grid on

    % Calculate implied volatilities from Lewis formula prices
    [impvols_Lewis] = plotVolSmile(pricesLewis, moneyness, r, T, S0);

    % Calculate implied volatilities from Exact Decomposition prices
    [impvols_ED] = plotVolSmile(prices_ED, moneyness, r, T, S0);

    % Calculate implied volatilities from Fast General Monte Carlo prices
    [impvols_FM] = plotVolSmile(prices_FM, moneyness, r, T, S0);

    % Calculate implied volatilities from confidence interval bounds of Exact Decomposition
    [impvols_ED_low] = plotVolSmile(conf_int_ED(1,:), moneyness, r, T, S0);
    [impvols_ED_up] = plotVolSmile(conf_int_ED(2,:), moneyness, r, T, S0);

    % Calculate implied volatilities from confidence interval bounds of Fast General Monte Carlo
    [impvols_FM_low] = plotVolSmile(conf_int_FM(1,:), moneyness, r, T, S0);
    [impvols_FM_up] = plotVolSmile(conf_int_FM(2,:), moneyness, r, T, S0);

    legend('Prices with Lewis', 'Mid Derivative Prices ED', 'Mid Derivative Prices FM', ...
           'Lower Bound ED', 'Upper Bound ED', 'Lower Bound FM', 'Upper Bound FM')

end % function plotVolatilities
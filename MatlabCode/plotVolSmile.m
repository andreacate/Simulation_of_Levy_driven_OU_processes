function [impvols] = plotVolSmile(prices, moneyness, rate, T, S0)
% Plots of volatility smile
%
% INPUT:
% prices:       prices of option
% moneyness:    values of moneyness
% T:            time to maturity
% S0:           initial underlying condition

% OUTPUTS:
% impvols:  A vector containing the calculated implied volatilities.

    %% Calculate strikes from moneyness
    % calculation of the corresponding strike prices for the given moneyness values.
    strikes = exp(-moneyness) .* S0 .* exp(T * rate);

    %% Calculate implied volatilities with Black-Scholes model
    % Calculation of Black implied volatilities given prices
    impvols = blkimpv(S0 * exp(T * rate), strikes, rate, T, prices);

    %% Plot the implied volatility smile

    plot(strikes, impvols);

end % function plotVolSmile


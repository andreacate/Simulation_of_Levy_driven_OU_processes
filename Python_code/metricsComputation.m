function [normInf, RMSE, MAPE] = metricsComputation (pricesExact, prices)
% Computation of evaluation metrics
%
% INPUT:
% pricesExact:    prices with Lewis formula
% prices:         computed prices
%
% OUTPUT:
% normInf:        infinite norm
% RMSE:           root mean squared error
% MAPE:           mean absolute percentage error


    %% Maximum error (Norm Inf)
    normInf = norm(pricesExact - prices, 'Inf');

    %% Root Mean Squared Error
    mse = sum((pricesExact-prices).^2)/length(pricesExact);
    RMSE = sqrt(mse);

    %% Mean Absolute Percentage Error
    MAPE = 100 * 1/length(pricesExact) * sum( abs( (pricesExact-prices)./pricesExact ) );

end % function metricsComputation
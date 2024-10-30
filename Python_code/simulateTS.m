function [S, S_mean] = simulateTS(alpha, beta, c, Nsim)
% Simulate a Tempered Stable random variable, as described in:
% Y. Qu, A. Dassios & H. Zhao, 2021,
%      "Exact Simulation of Ornstein-Uhlenbeck Tempered Stable Processes"
%       in Appendix D.1,
% i.e. with a method based on the "SSR method".
% The acceptance/rejection algorithm is meant to simulate only the
% following part of the characteristic function:
%       exp(   c * (beta-1i*u).^ alpha    * gamma(-alpha) ...
%            - c *  beta       ^ alpha    * gamma(-alpha) );
% The obtained distribution is a fully positive-supported distribution,
% i.e., de facto, a TS subordinator.
%
% INPUT
%   alpha, beta, c: parameters of TS Lévy measure
%   Nsim: number of samples to generate
%
% OUTPUT
%   S:        vector of samples from the TS distribution
%   S_mean:   mean of the distribution
%

    %% Acceptance Rejection Simulation

    % Initialize the flag array and the output array S
    flag = zeros(Nsim, 1);
    S = zeros(Nsim, 1);

    while sum(flag) < Nsim
        idx = find(flag == 0); % Indices of the remaining samples to be generated
        % [1] generate a stable random variable via Zolotarev
        Us = (rand(size(idx)) - 0.5) * pi;
        Es = -log(rand(size(Us))); % generating exponential r.v. (exprnd is slow!)

        S_temp = (-c * gamma(-alpha))^(1/alpha) * ...
            sin(alpha * Us + 0.5*pi*alpha) ./ cos(Us).^(1/alpha) .* ...
            (cos((1-alpha)*Us - 0.5*pi*alpha) ./ Es).^((1-alpha)/alpha);

        % [2] generate a uniform U
        U = rand(size(Us));

        % [3] A/R condition
        accept = (U <= exp(-beta * S_temp));

        % Update the flag and S arrays
        flag(idx(accept)) = 1;
        S(idx(accept)) = S_temp(accept);
    end

    %% Corrections

    % add "drift" part of the characteristic function:
    % the one associated with the third term of the Lévy-Khintchine
    S_mean = c * beta^(alpha-1) * gamma(-alpha+1);

end
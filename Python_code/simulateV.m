function v_vector = simulateV(alpha, a, nSim)
% Simulate the random variable V, introduced in:
% Sabino P. and Cufaro Petroni N., 2022, 
%   "Fast simulation of tempered stable Ornstein–Uhlenbeck processes"
% and used in the simulation of finite variation OU-TS and infinite
% variation OU-NTS.
%
% INPUT
%   modelParams: struct containing specifics for the selected process
%   a:           exp(-b*dt)
%   nSim:        number of simulations
%
% OUTPUT
%   v_vector:    vector of samples from the variable "V"
%

%% Extracting parameters

% driver parameters
% alpha = modelParams.alpha;

%% simulation

% number of intervals in which to divide [0,1]. 
L = 100; % as in the paper: Sabino P. & Cufaro-Petroni N., 2022

% define the pdf W, the modification of the original pdf V which is 
% convex and increasing
pdf_W = @(x) -a^alpha * log(a^alpha) / ...
    (1 - a^alpha + a^alpha * log(a^alpha)) * ...
    (a^(-alpha*x) - 1);

% generate a grid with L+1 ticks
w_grid = linspace(0, 1, L+1);

% and compute the corresponding values for the PDF
fw_grid = 0 * w_grid;
for ii = 1:length(w_grid)
    fw_grid(ii) = pdf_W( w_grid(ii) );
end

%% define the new rescaled quantities

% ql corresponds to the integral of the interval [w_l, w_{l+1}] (trapezio)
ql_grid = 0.5 * ( fw_grid(1:end-1) + fw_grid(2:end) ) / L;

% GL is the total area below the piecewise linear function
GL = sum(ql_grid);

% pl is the normalized ql, so that we create a new probability measure
pl_grid = ql_grid / GL;

% we compute the cdf of pl_grid, to simulate a discrete rv. We add 0 in
% front of the vector. The length of cumul_pl_grid is L.
cumul_pl_grid = cumsum(pl_grid);

%% A/R simulation

% initialization
v_vector = zeros(nSim,1);

for ii = 1:length(v_vector)

    % flag per acceptance/rejection
    accettato = false;

    while ~accettato

        % [1] we need to draw a rv s, i.e. to identify a small interval I_s
        s = find(cumul_pl_grid > rand(1), 1, 'first');
        % "s" defines the selected small interval. s=1 is the first 
        % interval, that corresponds to [w_grid(1), w_grid(2)]

        % [2] we need to simulate the rv g_s, corresponding to interval I_s
        coeff_ang = (fw_grid(s+1) - fw_grid(s)) / (w_grid(s+1) - w_grid(s));

        % and, as the pdf of g_s is obtained via linear interpolation,
        % I manually invert the CDF for the simulation
        y = w_grid(s) - fw_grid(s) / coeff_ang + ...
            sqrt(fw_grid(s)^2 / coeff_ang^2 + 2 / coeff_ang * ql_grid(s) * rand(1));

        % [3] I simulate the uniform rv "u", for the A/R method
        u = rand(1);

        % [4] A/R condition
        fw_y = pdf_W(y);
        gl_y = (coeff_ang * (y-w_grid(s)) + fw_grid(s));

        % [5] exit clause
        accettato = (u <= (fw_y/gl_y));

    end

    % I use the inverse transformation to obtain the original rv V
    % NB: there is a typo in the paper
    v_vector(ii) = a^(-y);

end


end
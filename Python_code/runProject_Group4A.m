% Final Project
% Simulation of Lévy-driven OU processes
% AY 2023-2024
% Group 4a - Catelli Andrea, Marchetto Erica, Urso Giovanni

% [1] Baviera & Manzoni, 2024, "Fast and General Simulation of L ́evy-driven 
% OU processes for Energy Derivatives"

% [2] Longstaff & Schwartz, 2001, "Valuing American Options by Simulation: 
% A Simple Least-Squares Approach"

% [3] Sabino, 2022, "Pricing Energy Derivatives in Markets Driven 
% by Tempered Stable and CGMY Processes of Ornstein–Uhlenbeck Type"

% [4] Sabino & Cufaro Petroni, 2022, "Fast simulation of tempered stable orn-
% stein–uhlenbeck processes"


clear all; close all; clc;
format long;


%% Data

% Model Parameters
x0 = 0;
b = 0.1;
beta_p = 2.5;
beta_n = 3.5;
c_p = 0.5;
c_n = 1;
gamma_c = 0;

% Flags
flag_OU_TS_FA = 1;
flag_TS_OU = 2;
flag_OU_TS_FV = 3;

flag_ED = 1;
flag_FGMC = 2;


%% Simulations parameters

T = 1/12;
Nsim = 1e7;
M = 1;


%% 1) TS-OU and OU-TS: exact simulation

%% OU-TS Finite Activity Exact

% Alpha definition
alpha = -1;

% Simulations
tic
[X_OU_TS_FA_ED, theorCumulantsT, simCumulantsT, theorCumulants_dt, simCumulants_dt, logFwd_OU_TS_FA_ED] = ...
                sim_OU_TS_FinAct_ED(x0, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c, T, M, Nsim);
toc

% Theoretical and simulated cumulants 
printCumulants (theorCumulantsT, simCumulantsT, theorCumulants_dt, simCumulants_dt, T, M)

% Check of martingality of the process
fprintf("Check of martingality: %d \n", mean(exp(logFwd_OU_TS_FA_ED)))


%% Plot of some simulations OU-TS Finite Activity

plotSimulation (x0, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c, 1, 1000, 50, 1)


%% OU-TS Finite Variation Exact

% Simulations
alpha = 0.5;

tic
[X_OU_TS_FV_ED, theorCumulantsT, simCumulantsT, theorCumulants_dt, simCumulants_dt, logFwd_OU_TS_FV_ED] = ...
                sim_OU_TS_FinVar_ED(x0, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c, T, M, Nsim);
toc

% Theoretical and simulated cumulants 
printCumulants (theorCumulantsT, simCumulantsT, theorCumulants_dt, simCumulants_dt, T, M)

% Check of martingality of the process
fprintf("Check of martingality: %d\n", mean(exp(logFwd_OU_TS_FV_ED)))


%% Plot of some simulations OU-TS Finite Variation
alpha=0.5;
plotSimulation (x0, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c, 1, 1000, 10, 3)
grid on

%% TS-OU Finite Variation Exact

% Simulations
alpha = 0.5;

tic
[X_TS_OU_FV_ED, theorCumulantsT, simCumulantsT, theorCumulants_dt, simCumulants_dt, logFwd_TS_OU_FV_ED] = ...
                sim_TS_OU_FinVar_ED(x0, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c, T, M, Nsim);
toc

% Theoretical and simulated cumulants 
printCumulants (theorCumulantsT, simCumulantsT, theorCumulants_dt, simCumulants_dt, T, M)

% Check of martingality of the process
fprintf("Check of martingality: %d \n", mean(exp(logFwd_TS_OU_FV_ED)))


%% Plots of some simulations TS-OU Finite Variation

plotSimulation (x0, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c, 1, 1000, 10, 2)


%% 2a) Fast MC Simulation

%% OU-TS Finite Variation FGMC

% Simulations
M_fft = 16;
scale = 1;
alpha = 0.5;

tic
[X_OU_TS_FV_FM, theorCumulantsT, simCumulantsT, theorCumulants_dt, simCumulants_dt, logFwd_OU_TS_FV_FM] = ...
                                             sim_OU_TS_FinVar_FGMC (x0, b, alpha, ...
                                             beta_p, beta_n, c_p, c_n, gamma_c, T, Nsim, M, M_fft, scale);
toc
    
% Theoretical and simulated cumulants 
printCumulants (theorCumulantsT, simCumulantsT, theorCumulants_dt, simCumulants_dt, T, M)

% Check of martingality of the process
fprintf("Check of martingality: %d \n", mean(exp(logFwd_OU_TS_FV_FM)))


%% TS-OU Finite Variation FGMC

% Simulations 
M_fft = 24;
scale = 1;
alpha = 0.5;

tic
[X_TS_OU_FV_FM, theorCumulantsT, simCumulantsT, theorCumulants_dt, simCumulants_dt, logFwd_TS_OU_FV_FM] = ...
                                 sim_TS_OU_FinVar_FGMC (x0, b, alpha, ...
                                 beta_p, beta_n, c_p, c_n, gamma_c, T, Nsim, M, M_fft, scale);     
toc

% Theoretical and simulated cumulants 
printCumulants (theorCumulantsT, simCumulantsT, theorCumulants_dt, simCumulants_dt, T, M)

% Check of martingality of the process
fprintf("Check of martingality: %d \n", mean(exp(logFwd_TS_OU_FV_FM)))


%% 2b) Fast MC Simulation for Finite Activity

%% OU-TS Finite Activity FGMC

% Simulations parameters
M_fft = 16;
scale = 1;
alpha = -1;

tic
[X_OU_TS_FA_FM, theorCumulantsT, simCumulantsT, theorCumulants_dt, simCumulants_dt, logFwd_OU_TS_FA_FM] = ...
                                 sim_OU_TS_FinAct_FGMC (x0, b, alpha, ...
                                 beta_p, beta_n, c_p, c_n, gamma_c, T, Nsim, M, M_fft, scale);
toc       

% Theoretical and simulated cumulants 
printCumulants (theorCumulantsT, simCumulantsT, theorCumulants_dt, simCumulants_dt, T, M)

% Check of martingality of the process
fprintf("Check of martingality: %d \n", mean(exp(logFwd_OU_TS_FA_FM)))


%% 3a) Energy derivative pricing - EU CALL

%% OU-TS Finite Variation

% Quantities of interest
alpha = 0.5;
S0 = 1; 
r = 0;
scale = 1;
T = 1/12;
M_fft = 16;
moneyness = linspace(-0.2,0.2,10)*sqrt(T);

% Prices with Lewis formula
tic
pricesLewis = priceEuropeanLewis_FFT(S0, b, alpha, beta_p, beta_n, c_p, c_n, gamma_c, T, moneyness, r, flag_OU_TS_FV, M_fft);
toc

% Prices exact decomposition
tic
[prices_ED, SD_ED, conf_int_ED] = priceEuropean (x0, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c, T, M, M_fft, Nsim, r, S0, moneyness, scale, flag_OU_TS_FV, flag_ED);
toc

% Prices with Fast General Monte Carlo
tic
[prices_FM, SD_FM, conf_int_FM] = priceEuropean (x0, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c, T, M, M_fft, Nsim, r, S0, moneyness, scale, flag_OU_TS_FV, flag_FGMC);
toc

% Prices in percentage 
pricesLewis_100 = pricesLewis*100;
prices_ED_100 = prices_ED*100;
prices_FM_100 = prices_FM*100;


%% Error computation

[normInf_ED, RMSE_ED, MAPE_ED] = metricsComputation (pricesLewis, prices_ED);
fprintf("Metrics value Exact decomposition \n")
fprintf("Norm Infinitive : %d \n", normInf_ED)
fprintf("Root Mean Squared Error : %d \n", RMSE_ED)
fprintf("Mean Absolute Percentage Error : %d \n", MAPE_ED)
fprintf("Mean Standard Deviation: %d \n", mean(SD_ED))

[normInf_FM, RMSE_FM, MAPE_FM] = metricsComputation (pricesLewis, prices_FM);
fprintf("Metrics value FGMC \n")
fprintf("Norm Infinitive : %d \n", normInf_FM)
fprintf("Root Mean Squared Error : %d \n", RMSE_FM)
fprintf("Mean Absolute Percentage Error : %d \n", MAPE_FM)
fprintf("Mean Standard Deviation: %d \n", mean(SD_FM))


%% Plot
    
plotPrices(moneyness, pricesLewis, prices_ED, conf_int_ED, prices_FM, conf_int_FM)


%% Implied Volatility

plotVolatilities (pricesLewis, prices_ED, prices_FM, conf_int_ED, conf_int_FM, moneyness, r, T, S0);

 
%% OU-TS Finite Activity 

% Quantities of interest
alpha = -1;
M_fft = 16;
S0 = 1; 
T = 1/12;
scale = 1;
r = 0;
moneyness = linspace(-0.2,0.2,10) * sqrt(T);

% Prices with Lewis formula
tic
pricesLewis = priceEuropeanLewis_FFT(S0, b, alpha, beta_p, beta_n, c_p, c_n, gamma_c, T, moneyness, r, flag_OU_TS_FA, M_fft);
toc

% Prices exact decomposition
tic
[prices_ED, SD_ED, conf_int_ED] = priceEuropean (x0, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c, T, M, M_fft, Nsim, r, S0, moneyness, scale, flag_OU_TS_FA, flag_ED);
toc

% Prices with Fast General Monte Carlo
tic
[prices_FM, SD_FM, conf_int_FM] = priceEuropean (x0, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c, T, M, M_fft, Nsim, r, S0, moneyness, scale, flag_OU_TS_FA, flag_FGMC);
toc

% Prices in percentage
pricesLewis_100 = pricesLewis*100;
prices_ED_100 = prices_ED*100;
prices_FM_100 = prices_FM*100;


%% Error computation

[normInf_ED, RMSE_ED, MAPE_ED] = metricsComputation (pricesLewis, prices_ED);
fprintf("Metrics value Exact decomposition \n")
fprintf("Norm Infinitive : %d \n", normInf_ED)
fprintf("Root Mean Squared Error : %d \n", RMSE_ED)
fprintf("Mean Absolute Percentage Error : %d \n", MAPE_ED)

[normInf_FM, RMSE_FM, MAPE_FM] = metricsComputation (pricesLewis, prices_FM);
fprintf("Metrics value FGMC \n")
fprintf("Norm Infinitive : %d \n", normInf_FM)
fprintf("Root Mean Squared Error : %d \n", RMSE_FM)
fprintf("Mean Absolute Percentage Error : %d \n", MAPE_FM)

%% Plot

plotPrices(moneyness, pricesLewis, prices_ED, conf_int_ED, prices_FM, conf_int_FM)


%% Implied Volatility

plotVolatilities (pricesLewis, prices_ED, prices_FM, conf_int_ED, conf_int_FM, moneyness, r, T, S0);

 
%% TS-OU Finite Variation 

% Quantities of interest
alpha = 0.5;
M_fft = 24;
S0 = 1; 
r = 0;
T = 1/12;
moneyness = linspace(-0.2,0.2,10)*sqrt(T);
scale = 1;

% Prices with Lewis formula using M_fft = 16
tic
pricesLewis = priceEuropeanLewis_FFT(S0, b, alpha, beta_p, beta_n, c_p, c_n, gamma_c, T, moneyness, r, flag_TS_OU, 16);
toc

% Prices exact decomposition
tic
[prices_ED, SD_ED, conf_int_ED] = priceEuropean (x0, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c, T, M, M_fft, Nsim, r, S0, moneyness, scale, flag_TS_OU, flag_ED);
toc

% Prices with Fast General Monte Carlo
tic
[prices_FM, SD_FM, conf_int_FM] = priceEuropean (x0, alpha, b, beta_p, beta_n, c_p, c_n, gamma_c, T, M, M_fft, Nsim, r, S0, moneyness, scale, flag_TS_OU, flag_FGMC);
toc


%% Error computation

[normInf_ED, RMSE_ED, MAPE_ED] = metricsComputation (pricesLewis, prices_ED);
fprintf("Metrics value Exact decomposition \n")
fprintf("Norm Infinitive : %d \n", normInf_ED)
fprintf("Root Mean Squared Error : %d \n", RMSE_ED)
fprintf("Mean Absolute Percentage Error : %d \n", MAPE_ED)

[normInf_FM, RMSE_FM, MAPE_FM] = metricsComputation (pricesLewis, prices_FM);
fprintf("Metrics value FGMC \n")
fprintf("Norm Infinitive : %d \n", normInf_FM)
fprintf("Root Mean Squared Error : %d \n", RMSE_FM)
fprintf("Mean Absolute Percentage Error : %d \n", MAPE_FM)

%% Plot

plotPrices(moneyness, pricesLewis, prices_ED, conf_int_ED, prices_FM, conf_int_FM)


%% Implied Volatility

plotVolatilities (pricesLewis, prices_ED, prices_FM, conf_int_ED, conf_int_FM, moneyness, r, T, S0);


%% 3b) Energy derivative pricing - AMERICAN CALL

% Quantities of interest
T = 1;
M = 52;
r = 0;
moneyness = linspace(-0.2,0.2,5)*sqrt(T);
K =  exp(r*T) * exp(-moneyness);
S0 = 1;
Nsim = 1e6;


%% OU - TS Finite Variation

alpha = 0.5;
M_fft = 18;
scale = 1;

% Exact Decomposition
tic
[price_ED, conf_int_ED] = priceAmerican (x0, alpha, b, beta_p, beta_n, ...
    c_p, c_n, gamma_c, T, M, M_fft, Nsim, K, r, S0, scale, flag_OU_TS_FV, flag_ED);
toc

% Fast General MC
tic
[price_FM, conf_int_FM] = priceAmerican (x0, alpha, b, beta_p, beta_n, ...
    c_p, c_n, gamma_c, T, M, M_fft, Nsim, K, r, S0, scale, flag_OU_TS_FV, flag_FGMC);
toc 

% Lewis formula
tic
pricesLewis = priceEuropeanLewis_FFT(S0, b, alpha, beta_p, beta_n, ...
    c_p, c_n, gamma_c, T, moneyness, r, flag_OU_TS_FV, M_fft);
toc


%% Plot

plotPrices(moneyness, pricesLewis, price_ED, conf_int_ED, price_FM, conf_int_FM)


%% 4) Gaussian OU

%% Process simulation Exact
x0 = 0;
b = 0.1;
sigma = 0.2;
M = 1;
Nsim = 1e7;
T = 1/12;

%Exact Decomposition Gaussian
tic
[X, theorCumulantsT, simCumulantsT, theorCumulants_dt, simCumulants_dt, logFwd_G] = sim_Gaussian_OU_ED(x0, b, sigma, T, M, Nsim);
toc

% Theoretical and simulated cumulants 
printCumulants (theorCumulantsT, simCumulantsT, theorCumulants_dt, simCumulants_dt, T, M)

% Check of martingality of the process
fprintf("Check of martingality: %d \n", mean(exp(logFwd_G)))


%% Pricing with Gaussian OU European

Nsim = 1e7;
T = 1/12;
M = 1;
moneyness = linspace(-0.2,0.2,10) * sqrt(T);
r = 0;
b = 0.1;
sigma = 0.2;
S0 = 1;
M_fft = 16;

% Prices with Lewis formula
tic
pricesLewis = priceEuropeanLewis_FFT_gaussian(S0, b, sigma,T, moneyness, r,  M_fft);
toc

% Computed prices
tic
[prices_G, std_dev, conf_int_G] = priceEuropeanGaussian(S0, b, sigma, T, r, moneyness, Nsim);
toc

%% Plots

figure()
plot(moneyness, prices_G)
hold on
grid on
plot(moneyness, conf_int_G(1,:))
plot(moneyness, conf_int_G(2,:))
plot(moneyness, pricesLewis,'*-')
legend('Prices Gaussian OU', 'Lower Bound', 'Upper Bound', 'Prices with Lewis','FontSize',15)
title('Gaussian OU European Prices')


%% Pricing with Gaussian OU American

Nsim = 1e6;
M = 52;
T = 1;
moneyness = linspace(-0.2,0.2,5) * sqrt(T);
r = 0;
b = 0.1;
sigma = 0.2;
S0 = 1;
M_fft = 16;

% Prices American
tic
[prices_G_A, conf_int_G_A] = priceAmericanGaussian(S0, b, sigma, T, r, moneyness, M, Nsim);
toc

% Prices with Lewis formula the european
tic
pricesLewis = priceEuropeanLewis_FFT_gaussian(S0, b, sigma,T, moneyness, r,  M_fft);
toc

%% Plots

figure()
plot(moneyness, prices_G_A)
hold on
grid on
plot(moneyness, conf_int_G_A(1,:))
plot(moneyness, conf_int_G_A(2,:))
plot(moneyness, pricesLewis,'*-')
legend('Prices Gaussian OU', 'Lower Bound', 'Upper Bound', 'Prices with Lewis','FontSize',15)
title('Gaussian OU American Prices')


%% 5) TS-OU Finite Variation Calibration

% Quantities of interest
T = 1/12;
M = 1;
alpha = 0.5;
M_fft = 16;
S0 = 1; 
r = 0;
moneyness = linspace(-0.2,0.2,100) * sqrt(T);

% Prices with Lewis formula
pricesLewis_OU_TS = priceEuropeanLewis_FFT(S0, b, alpha, beta_p, beta_n, c_p, c_n, gamma_c, T, moneyness, r, flag_OU_TS_FV, M_fft);

% Calibration of the optimal parameters
[b_opt, alpha_opt, beta_p_opt, beta_n_opt, c_p_opt, c_n_opt, gamma_c_opt] = ...
                   calibration(S0, T, moneyness, r, flag_TS_OU, M_fft, pricesLewis_OU_TS);

               
%% Prices of TS-OU with Lewis formula with the new optimal parameters

pricesLewis_calibrated = priceEuropeanLewis_FFT(S0, b_opt, alpha_opt, beta_p_opt, beta_n_opt, c_p_opt, c_n_opt, gamma_c_opt, T, moneyness, r, flag_TS_OU, M_fft);

figure()
plot(moneyness, pricesLewis_OU_TS)
hold on
grid on
plot(moneyness, pricesLewis_calibrated)
legend('Prices Lewis OU-TS', 'Prices Lewis TS-OU Calibrated','FontSize', 15)
title('Prices Lewis OU-TS vs TS-OU')

%% Implied Volatility with the new optimal parameters

figure()
impvols_calibrated = plotVolSmile(pricesLewis_calibrated, moneyness, r, T, S0);
hold on
grid on
impvols_OU_TS = plotVolSmile(pricesLewis_OU_TS, moneyness, r, T, S0);
legend('Implied Volatility TS-OU Calibrated', 'Implied Volatility OU-TS','FontSize', 15)
title('Volatility Smiles')



%% Additional analysis: American Call with the other processes

% %% OU - TS Finite Activity
% 
% alpha = -1;
% M_fft = 16;
% scale = 1;
% 
% % Exact Decomposition
% tic
% [price_ED, conf_int_ED] = priceAmerican (x0, alpha, b, beta_p, beta_n, ...
%     c_p, c_n, gamma_c, T, M, M_fft, Nsim, K, r, S0, scale, flag_OU_TS_FA, flag_ED);
% toc
% 
% % Fast General MC
% tic
% [price_FM, conf_int_FM] = priceAmerican (x0, alpha, b, beta_p, beta_n, ...
%     c_p, c_n, gamma_c, T, M, M_fft, Nsim, K, r, S0, scale, flag_OU_TS_FA, flag_FGMC);
% toc
% 
% % Lewis formula
% pricesLewis = priceEuropeanLewis_FFT(S0, b, alpha, beta_p, beta_n, ...
%     c_p, c_n, gamma_c, T, moneyness, r, flag_OU_TS_FV, M_fft);
% 
% 
% %% Plot
% 
% plotPrices(moneyness, pricesLewis, price_ED, conf_int_ED, price_FM, conf_int_FM)
% 
% 
% %% TS - OU Finite Variation
% 
% alpha = 0.5;
% M_fft = 16;
% scale = 0.5;
% 
% % Exact Decomposition
% tic
% [price_ED, conf_int_ED] = priceAmerican (x0, alpha, b, beta_p, beta_n, ...
%     c_p, c_n, gamma_c, T, M, M_fft, Nsim, K, r, S0, scale, flag_TS_OU, flag_ED);
% toc
% 
% % Fast General MC
% tic
% [price_FM, conf_int_FM] = priceAmerican (x0, alpha, b, beta_p, beta_n,...
%     c_p, c_n, gamma_c, T, M, M_fft, Nsim, K, r, S0, scale, flag_TS_OU, flag_FGMC);
% toc
% 
% % Lewis formula
% pricesLewis = priceEuropeanLewis_FFT(S0, b, alpha, beta_p, beta_n,...
%     c_p, c_n, gamma_c, T, moneyness, r, flag_OU_TS_FV, M_fft);
% 
% 
% %% Plot
% 
% plotPrices(moneyness, pricesLewis, price_ED, conf_int_ED, price_FM, conf_int_FM)

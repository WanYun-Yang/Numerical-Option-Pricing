%% Final project for scientific computing _ Black-Scholes equation

%  We want to calculate this equation in the following three ways to 
%  compare the performance and accuracy.

% Chien-Hao Wu, Wan-Yun Yang, Ting-Wei Su, Dec. 2017. Last changed, 23/01/2018

%% Inputs
% Smax:   the maximum possible stock price
% S0:     the initial stock price
% K:      the strike price
% r:      the risk-free interest rate
% sigma:  volatility of the market
% T:      time to maturity
% n:      the number of time intervals
% m:      the number of spaing for the stock price
% type:   Call or put, American or European

Smax = 100; 
printmax = Smax;
S0 = 20;
K = 21;
r = 0.1;
sigma = 0.3;
T = 4/12;
n = 600;
m = 600;

%% Solutions
% Use forward difference to solve(discrete directly)
EuCall_FD = FD_OptionPricing(Smax, T, K, r, m, n, S0, sigma, 'EC'); 
EuPut_FD = FD_OptionPricing(Smax, T, K, r, m, n, S0, sigma, 'EP');
AmCall_FD = FD_OptionPricing(Smax, T, K, r, m, n, S0, sigma, 'AC'); 
AmPut_FD = FD_OptionPricing(Smax, T, K, r, m, n, S0, sigma, 'AP');

% Use Crank-Nicolson method to solve(discrete directly)
EuCall_CN = CN_OptionPricing(Smax, T, K, r, m, n, S0, sigma, 'EC', printmax);
EuPut_CN = CN_OptionPricing(Smax, T, K, r, m, n, S0, sigma, 'EP', printmax);
AmCall_CN = CN_OptionPricing(Smax, T, K, r, m, n, S0, sigma, 'AC', printmax); 
AmPut_CN = CN_OptionPricing(Smax, T, K, r, m, n, S0, sigma, 'AP', printmax);

% Use Crank-Nicolson method to solve(transform it into heat equation)
EuCall_HT = HT_OptionPricing(Smax, T, K, r, m, n, S0, sigma, 'EC');
EuPut_HT = HT_OptionPricing(Smax, T, K, r, m, n, S0, sigma, 'EP');
AmCall_HT = HT_OptionPricing(Smax, T, K, r, m, n, S0, sigma, 'AC'); 
AmPut_HT = HT_OptionPricing(Smax, T, K, r, m, n, S0, sigma, 'AP');

% Use Binomial Model Derivation
EuCall_BT = BT_OptionPricing(T, K, r, S0, sigma, 3600, 'EC');
EuPut_BT = BT_OptionPricing(T, K, r, S0, sigma, 3600, 'EP');
AmCall_BT = BT_OptionPricing(T, K, r, S0, sigma, 3600, 'AC'); 
AmPut_BT = BT_OptionPricing(T, K, r, S0, sigma, 3600, 'AP');

% use Build_in function to solve
[EuCall_BU, EuPut_BU] = blsprice(S0, K, r, T, sigma);  

% Show the result
Type = {'EuCall';'EuPut';'AmCall';'AmPut'};
Directly_FD = [EuCall_FD; EuPut_FD; AmCall_FD; AmPut_FD];
Directly_CN = [EuCall_CN; EuPut_CN; AmCall_CN; AmPut_CN];
Heat_CN = [EuCall_HT; EuPut_HT; AmCall_HT; AmPut_HT];
BinomialTree = [EuCall_BT; EuPut_BT ;AmCall_BT; AmPut_BT];
Table = table(Directly_FD,Directly_CN, Heat_CN, BinomialTree, 'RowNames', Type)






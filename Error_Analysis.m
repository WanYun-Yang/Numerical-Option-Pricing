%% Error Analysis
%% Settings
Smax = 100; 
printmax = Smax;
K = 21;
r = 0.1;
sigma = 0.3;
T = 4/12;
n = 50;
m = 50;
% Exactsol = @(s) normcdf((log(s/K) + (r + sigma^2 /2)*T)/sigma/sqrt(T)).*s ...
%                 - normcdf((log(s/K) + (r + sigma^2 /2)*T)/sigma/sqrt(T) - sigma*sqrt(T))*K*exp(-r*(T));
% S = Smax/m:Smax/m:Smax - Smax/m;
 
%% Error performance for different Smax
figure(1)
subplot(1, 2, 1)
len = 5;
error = zeros(2, len);
CNS1 = CN_E(Smax*2^len, T, K, r, m*2^len, n, sigma, 'EC');
CNS1 = CNS1(1:m-1);
HTS1 = HT_E(Smax*2^len, T, K, r, m*2^len, n, sigma, 'EC');
HTS1 = HTS1(1:m-1);
for i = 1:len
    CN1 = CN_E(Smax, T, K, r, m*2^(i - 1), n, sigma, 'EC');
    CN1 = CN1(1:m-1);
    error(1, i) = max(abs(CN1 -CNS1));
    HT1 = HT_E(Smax, T, K, r, m*2^(i - 1), n, sigma, 'EC');
    HT1 = HT1(1:m-1);
    error(2, i) = max(abs(HT1 - HTS1));
    Smax = Smax*2;
end

error = -log(error);
hold on

plot(1:len, error(1, :))
plot(1:len, error(2, :))
title('Error Performance for Different Smax')
xlabel('Smax*(2^i)')
ylabel('max Error(-log)')
legend('Direct-CN','Heat-CN')
Smax = 100;

%% Error performance for different delta T
n = 700;
m = 50;
len = 10;
error2 = zeros(3, len-1);
CNS = CN_E(Smax, T, K, r, m, n*2^len, sigma, 'EC');
HTS = HT_E(Smax, T, K, r, m, n*2^len, sigma, 'EC');
FDS = FD_E(Smax, T, K, r, m, n*2^len, sigma, 'EC');

for i = 1:len
    error2(1, i) = max(abs(CN_E(Smax, T, K, r, m, n, sigma, 'EC') - CNS));
    error2(2, i) = max(abs(HT_E(Smax, T, K, r, m, n, sigma, 'EC') - HTS));
    error2(3, i) = max(abs(FD_E(Smax, T, K, r, m, n, sigma, 'EC') - FDS));
    n = n*2;
end

error2 = -log(error2);
subplot(1, 2, 2)
hold on
plot(1:len, error2(1, :))
plot(1:len, error2(2, :))
plot(1:len, error2(3, :))
title('Error Performance for Different delta t')
xlabel('n*(2^i)')
ylabel('max Error(-log)')
legend('D-CN','Heat-CN', 'FTCS')

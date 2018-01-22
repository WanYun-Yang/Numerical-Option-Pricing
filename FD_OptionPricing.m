function sol= FD_OptionPricing(Smax, T, K, r, m, n, S0, sigma, type)
    % FD_OptionPricing for pricing European call and put
    % ------------------------------------------------------------------------
    % Input:
    %   Smax    --- a sufficiently high stock price in the given time horizon
    %   T       --- the maturity time of the option expressed in years
    %   K       --- the strike/exercise price
    %   r       --- the risk-free interest rate expressed continuously
    %               compounded per annum
    %   m       --- spacing in the spatial directiom (stock price) S
    %   n       --- spacing in the temporal direction t
    %   S0      --- the intial/current stock price
    %   sigma   --- the stock price volatility expressed per annum
    %   type    --- string character indicating which type of options is being
    %               considered, 'EC' for European call and 'EP' for European put
    %-------------------------------------------------------------------------
    % Output:
    %   sol     --- the option price
    %-------------------------------------------------------------------------
    % Example:
    % To value an European call option with 
    % T=4 months, K=21, S0=20, r=0.1, sigma=0.3, yield=0
    % use the command FD_OptionPricing(100, 4/12, 21, 0.1, 150, 150, 20, 0.3, 'EC').
    % Here I use n=150 and m=150 for spacing in temporal and spatial direction
    %
    % Chien-Hao Wu, Wan-Yun Yang, Dec. 2017. Last changed, 31/12/2017
    % ------------------------------------------------------------------------
    %%  Construct initial conditions
   % Step size for temporal and spatial direction
   n = ceil(T*(r + (sigma*m)^2));
   deltaT = T/n;
   deltaS = Smax/m;
     %% Boundary condition
    f0 = zeros(m+1, n+1);
    [bT, bS0, bSinf] = boundary(Smax, T, K, r, m, n, type);
    f0(1,:) = bS0;
    f0(m+1,:) = bSinf;
    f0(:,n+1) = bT;
    %% Solve the linear system
    
    alpha = sigma^2* deltaT/2.* (1:m-1)'.^2;
    beta = r.* (1:m-1)'.* deltaT/2;
    S = 0:Smax/m:Smax;
    if type == 'AC'
        for i = n:-1:1
            f0(2:end-1, i) = (alpha + beta).* f0(3:end, i+1)...
                            +(alpha - beta).* f0(1:end-2, i+1)...
                            +(1- r* deltaT - 2*alpha).* f0(2:end-1, i+1);
            f0(2:end-1, i) = max(f0(2:end-1, i), S(2:end-1)' - K);
        end
    elseif type == 'AP'
        for i = n:-1:1
            f0(2:end-1,i) = (alpha + beta).* f0(3:end, i+1)...
                            +(alpha - beta).* f0(1:end-2, i+1)...
                            +(1- r* deltaT - 2*alpha).* f0(2:end-1, i+1);
            f0(2:end-1, i) = max(f0(2:end-1, i), K - S(2:end-1)');
        end
    else 
        for i = n:-1:1
            f0(2:end-1,i) = (alpha + beta).* f0(3:end, i+1)...
                            +(alpha - beta).* f0(1:end-2, i+1)...
                            +(1- r* deltaT - 2*alpha).* f0(2:end-1, i+1);
        end
    end
    %% Compute the result
     % Interpolation
     F = f0;
     z = S0/ deltaS + 1;
     sol = double(F(floor(z):floor(z) + 1, 1)'*[floor(z) + 1 - z, z - floor(z)]');
end
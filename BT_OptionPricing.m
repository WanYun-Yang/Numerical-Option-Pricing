function sol = BT_OptionPricing(T, K, r, S0, sigma, NumT, type)
    % Inputs
    % S0: the initial stock price
    % K: the strike price
    % r: the risk-free interest rate
    % q: the yield on the underlying asset
    % sigma: volatility of the market
    % T: time to maturity in years
    % NumT: number of step in one year (must be a multiple of 12)
    % OptionType: 'p' for put and 'c' for call
    % ExerciseType: 'a' for American and 'e' for European
    % Outputs
    % price: the price of the option
    %-----------------------------------------------------------------------
    % Example:
    % To value an European call option with 
    % T=4 months, K=21, S0=20, r=0.1, sigma=0.3, yield=0
    % use the command my_BSmodel(20,21,0.1,0,0.3,4/12,3600,'c','e')
    %-----------------------------------------------------------------------

    NT = NumT*T;                %total number of time interval
    dt = T./NT;
    u = exp(sigma.*sqrt(dt));      %the up movement
    d = 1./u;                      %the down movement
    a = exp(r*dt);
    p = (a-d)./(u-d);              %the probability of an up movement
    f = zeros(NT+1,NT+1);          %initialization of the values
    
    % Option prices at maturity
    if (type == 'EP') | (type == 'AP') 
        f(1:NT+1,NT+1) = max(K-S0.*(u.^(NT:-1:0)).*(d.^(0:NT)),0);
    else
        f(1:NT+1,NT+1) = max(S0.*(u.^(NT:-1:0)).*(d.^(0:NT))-K,0);
    end
    
    %backwards induction process
    switch type
        case 'EC'
            for j = NT:-1:1 
                f(1:j, j) = flipud(max(exp(-r.*dt).*(p.*flipud(f(1:j, j+1)) + (1-p).*flipud(f(2:j+1, j+1))), 0));
            end
        case 'EP'
            for j = NT:-1:1 
                f(1:j, j) = flipud(max(exp(-r.*dt).*(p.*flipud(f(1:j, j+1)) + (1-p).*flipud(f(2:j+1, j+1))), 0));
            end
        case 'AC'
            for j = NT:-1:1
                f(1:j, j) = max([flipud(max(exp(-r.*dt).*(p.*flipud(f(1:j, j+1)) + (1-p).*flipud(f(2:j+1, j+1))), 0)),...
                            max(S0.*(u.^(j-1:-1:0)).*(d.^(0:j-1)) - K, 0)']');
            end
        case 'AP'
            for j = NT:-1:1 
                f(1:j, j) = max([flipud(max(exp(-r.*dt).*(p.*flipud(f(1:j, j+1)) + (1-p).*flipud(f(2:j+1, j+1))), 0)),...
                            max(K - S0.*(u.^(j-1:-1:0)).*(d.^(0:j-1)), 0)']');
            end
    end
    sol = f(1,1);
end
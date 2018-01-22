function sol = HT_OptionPricing(Smax, T, K, r, m, n, S0, sigma, type)
    % Settings
    T0 = T;
    T = T*sigma^2 /2;
    Sm0 = Smax;
    Smax = log(Smax/K);
    Smin = -20;
    c = r/(sigma^2)*2;
    fun = @(x, tau)exp((c-1)/2*x + (c+1)^2/4*tau)/K;
    
    % Construct boundary conditions
    M = Smin:(Smax-Smin)/m:Smax;
    N = fliplr(0:T/n:T);
    switch type
        case 'EC'  % European Call
        b0 = max(K*exp(M) - K, 0);
        bS0 = zeros(1, n+1);
        bSinf = Sm0 - K*exp(-r*(N/(sigma^2)* 2));
        case 'EP'  % European Put
        b0 = max(K - K*exp(M), 0);
        bS0 = K*exp(-r*(N/(sigma^2)* 2));
        bSinf = zeros(1, n+1);
        case 'AC'  % American Call
        b0 = max(K*exp(M) - K, 0);
        bS0 = zeros(1, n+1);
        bSinf = Smax - K;
        case 'AP'  % American Put
        b0 = max(K - K*exp(M), 0);
        bS0 = K*ones(1, n+1);
        bSinf = zeros(1, n+1);
        case 'DC'
        B = 100;
        b0 = B*heaviside(K*exp(M) - K);
        bS0 = zeros(1, n+1);
        bSinf = B*exp(-r*(N/(sigma^2)* 2));
    end
    
    [X, Y] = meshgrid(N, M);
    U = zeros(m+1, n+1);
    U(1, :) = bS0;
    U(m+1, :) = bSinf;
    U(:, end) = b0;
    U = U.*fun(Y, X);
    beta = (T/n)/((Smax-Smin)/m)^2/2;
    U = fliplr(U);

    % Solve the system
    if type == 'AC'
        for i = 2:n+1
            So = beta*U(1:end-2, i-1)...
                + (1 - 2*beta)*U(2:end-1, i-1)...
                + beta*U(3:end, i-1);
            So(1) = So(1) + beta*U(1, i); 
            So(end) = So(end) + beta*U(end, i); 
            U(2:end-1, i) = Thomas(So, -beta*ones(1, length(So)),...
               (1 + 2*beta)*ones(1, length(So)), -beta*ones(1, length(So)));
           VV = max(K*exp(-(c-1)/2*M - (c+1)^2*(i-1)*T/n/4)'.* U(:, i),...
               max(K*exp(M) - K, 0)');
           U(:, i) = fun(M, (i-1)*T/n)'.*VV;
        end
    elseif type == 'AP'
        for i = 2:n+1
            So = beta*U(1:end-2, i-1)...
                + (1 - 2*beta)*U(2:end-1, i-1)...
                + beta*U(3:end, i-1);
            So(1) = So(1) + beta*U(1, i); 
            So(end) = So(end) + beta*U(end, i); 
            U(2:end-1, i) = Thomas(So, -beta*ones(1, length(So)),...
               (1 + 2*beta)*ones(1, length(So)), -beta*ones(1, length(So)));
           VV = max(K*exp(-(c-1)/2*M - (c+1)^2*(i-1)*T/n/4)'.* U(:, i),...
               max(K - K*exp(M), 0)');
           U(:, i) = fun(M, (i-1)*T/n)'.*VV;
        end
    else
        for i = 2:n+1
            So = beta*U(1:end-2, i-1)...
                + (1 - 2*beta)*U(2:end-1, i-1)...
                + beta*U(3:end, i-1);
            So(1) = So(1) + beta*U(1, i); 
            So(end) = So(end) + beta*U(end, i); 
            U(2:end-1, i) = Thomas(So, -beta*ones(1, length(So)),...
               (1 + 2*beta)*ones(1, length(So)), -beta*ones(1, length(So)));
        end
    end
    
    V = K*exp(-(c-1)/2*Y - (c+1)^2*X/4).* fliplr(U);
    MM = K*exp(M);
    
    % Interpolation
    Stock = K*exp(M);
    [~, ind] = min(abs(Stock - S0));
    zz = spline(Stock, V(:, 1));
    if S0 < Stock(ind)
        z = S0 - Stock(ind-1);
        coe = zz.coefs(ind-1, :);
        sol = z^3*coe(1) + z^2*coe(2) + z*coe(3) + coe(4);
    elseif S0 > Stock(ind)
        z = S0 - Stock(ind);
        coe = zz.coefs(ind, :);
        sol = z^3*coe(1) + z^2*coe(2) + z*coe(3) + coe(4);
    else
        sol = V(ind, 1);
    end
end
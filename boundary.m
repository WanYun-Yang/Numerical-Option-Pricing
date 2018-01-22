function [bT, bS0, bSinf] = boundary(Smax, T, K, r, m, n, type)
    M = 0:Smax/m:Smax;
    N = 0:T/n:T;
    switch type
        case 'EC'  % European Call
            bT = max(M - K, 0);
            bS0 = zeros(1, n+1);
            bSinf = Smax - K*exp(-r*(T - N));
        case 'EP'  % European Put
            bT = max(K - M, 0);
            bS0 = K*exp(-r*(T - N));
            bSinf = zeros(1, n+1);
        case 'AC'  % American Call
            bT = max(M - K, 0);
            bS0 = zeros(1, n+1);
            bSinf = (Smax - K)* ones(1, n+1);
        case 'AP'  % American Put
            bT = max(K - M, 0);
            bS0 = K*ones(1, n+1);
            bSinf = zeros(1, n+1);
        case 'DC'  % Digital Call/Binary Call)
            B = 100;
            bT = B*heaviside(M - K);
            bS0 = zeros(1, n+1);
            bSinf = B*exp(-r*(T - N));
        case 'DP'  % Digital Put/Binary Put)
            B = 100;
            bT = B*(1-heaviside(M - K));
            bS0 = B*exp(-r*(T - N));
            bSinf = zeros(1, n+1);
    end
end
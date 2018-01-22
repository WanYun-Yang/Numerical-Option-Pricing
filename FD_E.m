function F0 = FD_E(Smax, T, K, r, m, n, sigma, type)
    %%  Construct initial conditions
   % Step size for temporal and spatial direction
   deltaT = T/n;
    
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
    
    F0 = f0(2:end-1, 1); 
end
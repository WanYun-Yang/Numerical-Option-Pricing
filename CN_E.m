function F0 = CN_E(Smax, T, K, r, m, n, sigma, type)
    % Construct boundary conditions
    F = zeros(m+1, n+1);
    [bT, bS0, bSinf] = boundary(Smax, T, K, r, m, n, type);
    F(1, :) = bS0;
    F(m+1, :) = bSinf;
    F(:, n+1) = bT;
    
    S = 0:Smax/m:Smax;
    alpha = (sigma*S*m/Smax/2).^2;
    beta = r* S /4 /(Smax/m);
    
    % Solve the system
    if type == 'AC'
        for j = n:-1:1
            U = (alpha(2: end-1) - beta(2: end-1))'.*F(1:end-2, j+1)...
                + (n/T - 2*alpha(2:end-1) - r/2)'.*F(2:end-1, j+1)...
                + (alpha(2: end-1) + beta(2: end-1))'.*F(3:end, j+1);
            U(1) = U(1) + (alpha(2) - beta(2))*F(1, j); 
            U(end) = U(end) + (alpha(end-1) + beta(end-1))*F(end, j); 
            F(2:end-1, j) = Thomas(U, beta(2: end-1) - alpha(2: end-1),...
                (n/T + 2*alpha(2: end-1) + r/2), - alpha(2: end-1) - beta(2: end-1));
            F(2:end-1, j) = max(F(2:end-1, j), S(2:end-1)' - K);
        end
    elseif type == 'AP'
        for j = n:-1:1
            U = (alpha(2: end-1) - beta(2: end-1))'.*F(1:end-2, j+1)...
                + (n/T - 2*alpha(2:end-1) - r/2)'.*F(2:end-1, j+1)...
                + (alpha(2: end-1) + beta(2: end-1))'.*F(3:end, j+1);
            U(1) = U(1) + (alpha(2) - beta(2))*F(1, j); 
            U(end) = U(end) + (alpha(end-1) + beta(end-1))*F(end, j); 
            F(2:end-1, j) = Thomas(U, beta(2: end-1) - alpha(2: end-1),...
                (n/T + 2*alpha(2: end-1) + r/2), - alpha(2: end-1) - beta(2: end-1));
            F(2:end-1, j) = max(F(2:end-1, j), K - S(2:end-1)' );
        end
    else
        for j = n:-1:1
            U = (alpha(2: end-1) - beta(2: end-1))'.*F(1:end-2, j+1)...
                + (n/T - 2*alpha(2:end-1) - r/2)'.*F(2:end-1, j+1)...
                + (alpha(2: end-1) + beta(2: end-1))'.*F(3:end, j+1);
            U(1) = U(1) + (alpha(2) - beta(2))*F(1, j); 
            U(end) = U(end) + (alpha(end-1) + beta(end-1))*F(end, j); 
            F(2:end-1, j) = Thomas(U, beta(2: end-1) - alpha(2: end-1),...
                (n/T + 2*alpha(2: end-1) + r/2), - alpha(2: end-1) - beta(2: end-1));
        end
    end

    F0 = F(2:end-1, 1); 
end
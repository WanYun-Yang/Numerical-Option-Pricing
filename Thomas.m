function x = Thomas(U, a, b, c)
    num = length(U);
    c1 = zeros(num, 1);
    d1 = c1;
    x = c1;
    c1(1) = c(1)/b(1);
    d1(1) = U(1)/b(1);
    for i = 2:num
        c1(i) = c(i)/(b(i) - a(i)*c1(i-1));
        d1(i) = (U(i) - a(i)*d1(i-1))/(b(i) - a(i)*c1(i-1));
    end
    x(end) = d1(end);
    for j = num-1:-1:1
        x(j) = d1(j) - c1(j)*x(j+1);
    end
end
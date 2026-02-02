function d = dS2(tau)
%2 sum[(-1)^n q^(n+0.5)^2 * (2n+1)pi]
    q = exp(pi*1i*(tau));
    d= 0;
    for n = 0:20
        diff = 2*(-1)^n * q^((n+0.5)*(n+0.5))*(2*n+1)*pi;
        if abs(diff) < 1e-20
            break
        end
        d = d+ diff;
    end
end

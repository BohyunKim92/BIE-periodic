function v1 = JacobiTheta1(z,tau)
% return JacobiTheta1 function value given z, complex number and tau = 1+ib, the elliptic period
% JacobiTheta1(z):= 2 * [sum_{n=0}^{\infty} (-1)^n * q^[(n+0.5)^2] * sin[(2n+1)*pi*z] ]

q = exp(pi*1i*(tau)); % nome
v1 = zeros(size(z));
for n = 0:20
    diff = 2*(-1)^n * q^((n+0.5)*(n+0.5))*sinpi((2*n+1)*z);
    if abs(diff) < 1e-20
        break
    end
    v1 = v1+ diff;
end
end
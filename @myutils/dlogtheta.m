function val = dlogtheta(z,tau)
% return dlog(JacobiTheta1(z))/dz given z, complex number and tau = 1+ib, the elliptic period
% dlog(JacobiTheta1(z))/dz calculated by (dJacobiTheta1(z)/dz)/JacobiTheta1(z)

q = exp(pi*1i*(tau)); % nome 
dv1 = zeros(size(z));

% calculating d/dz [JacobiTheta1(z)]
for n = 0:40
    diff = 2*(-1)^n * q^((n+0.5)^2)*((2*n+1)*pi)*cospi((2*n+1)*z);
    if abs(diff) < 1e-20
        break
    end
    dv1 =dv1+ diff;
end

% dividing by JacobiTheta1(z)
val = dv1./myutils.JacobiTheta1(z,tau);
end
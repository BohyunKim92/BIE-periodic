function [Gx,Gy] = gradG(z,tau)
% return [Gx,Gy] given z, complex number and tau = 1+ib, the elliptic period
% Gx = dG/dx = -1/(2*pi) * Re[log(JacobiTheta1(z))z]
% Gy = dG/dy = 1/(2*pi) * Im[log(JacobiTheta1(z))z] + Im[z]/b;  
% G(z):= -1/(2*pi)*log(abs(jacobi_theta1(z)))+1/(2*b)*Im(z)^2 for a complex number z. 

dlogv = myutils.dlogtheta(z,tau); % log(JacobiTheta1(z))z
b = imag(tau);
Gx = -1/(2*pi) * real(dlogv);
Gy = 1/(2*pi) * imag(dlogv)+ imag(z)./b;
end
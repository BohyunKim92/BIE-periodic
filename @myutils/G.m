function Geval = G(z,tau)
% return G(z) given z, complex number, and tau = 1+ib, the elliptic periods.
% G(z):= -1/(2*pi)*log(abs(jacobi_theta1(z)))+1/(2*b)*Im(z)^2 for a complex number z. 
b = imag(tau);
v1 = myutils.JacobiTheta1(z,tau); 
G1 = -1/(2*pi) * log(abs(v1));
G2 = 1/(2*b) * (imag(z).^2);
Geval = G1+G2;
end
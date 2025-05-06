function s = shiftedradial(M, fs, a1, varargin)
% create a radialfunction centered at a1.
f = fs{1};
Z = @(s) a1+exp(2i*pi*s).*f(2*pi*s);   % note conversion from 0<s<1 to 0<t<2pi
fp = fs{2};
if numel(fs)<=2
  s = segment(M, {Z, @(s) 2*pi*(1i*Z(s) + exp(2i*pi*s).*fp(2*pi*s))}, 'p');
else
  Zp = @(s) 2*pi*(1i*Z(s) + exp(2i*pi*s).*fp(2*pi*s));
  fpp = fs{3};
  Zpp = @(s) 4*pi^2*((fpp(2*pi*s) + 2i*fp(2*pi*s)).*exp(2i*pi*s) - Z(s));
  s = segment(M, {Z, Zp, Zpp}, 'p', varargin{:});
end
% Notice that all this routine did was manipulate function handles, no evals.


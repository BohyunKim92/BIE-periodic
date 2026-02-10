function [tval] = ztot(z,c)
zshift = z-c;
tval = zeros(size(zshift));
%given z, spit out value b/t 0 to 1, that represent angle [0,2pi]
q1 = @(z) logical((real(z) >=0).*(imag(z) >=0)); % first quadrant
q2 = @(z) logical((real(z) <0)); % second & third
q3 = @(z) logical((real(z) >=0).*(imag(z) <=0)); % 4th
ind1 = q1(zshift); zshift1 = zshift(ind1);
ind2 = q2(zshift); zshift2 = zshift(ind2);
ind3 = q3(zshift); zshift3 = zshift(ind3);

tval(ind1) = atan(imag(zshift1)./real(zshift1))./(2*pi);
tval(ind2) = (atan(imag(zshift2)./real(zshift2))+pi)./(2*pi);
tval(ind3) = (atan(imag(zshift3)./real(zshift3))+2*pi)./(2*pi);
end
function s = mycircle(c, r)
%MYCIRCLE Create a circular segment centered at a given complex location.
%   s = MYCIRCLE(c, r) returns a segment object representing a circle
%   centered at complex number 'c' with radius 'r'.

% Create the circular segment using the full parameterization [center, radius, start angle, end angle]
s = segment([], [c, r, 0, 2*pi]);

end

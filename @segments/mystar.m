function s = mystar(c, r, osc)
%MYSTAR Create a star-shaped segment centered at a given complex location.
%   s = MYSTAR(c, r, osc) returns a segment object with a star-like shape.
%   The shape is defined by an oscillatory radial function with frequency 'osc',
%   base radius 'r', and center at the complex location 'c'.

% Define radial function and its first two derivatives
fs = {
    @(t) 1 + r * cos(osc * t);             % radial function
    @(t) -osc * r * sin(osc * t);          % first derivative
    @(t) -osc^2 * r * cos(osc * t)         % second derivative
};

% Generate the star-shaped segment using 25 points
s = segment.shiftedradial(25, fs, 0);

% Normalize and translate the segment
s.scale(r / (r + 1));  % Scale to ensure radius is approximately 'r'
s.translate(c);        % Move the shape to center 'c'

end

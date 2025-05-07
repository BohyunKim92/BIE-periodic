function s = mytref(c, r)
%MYTREF Create a trefoil-shaped segment centered at a given complex point.
%   s = MYTREF(c, r) returns a segment object representing a trefoil shape
%   centered at complex location c with base radius r.
%
%   The shape is generated using a radial function with a cosine modulation
%   to create the trefoil lobes.

% Define trefoil radial function and its first two derivatives
tref  = @(q)  r * (1 + 0.3 * cos(3 * q));       % radial profile
dtref = @(q) -r * (0.3 * 3 * sin(3 * q));       % first derivative
d2tref = @(q) -r * (0.3 * 3^2 * cos(3 * q));    % second derivative

% Create segment using 25-point discretization and custom radial functions
s = segment.radialfunc(25, {tref, dtref, d2tref});

% Translate the segment to center c
s.translate(c);

end

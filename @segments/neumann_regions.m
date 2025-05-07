function [ss] = neumann_regions()
%NEUMANN_REGIONS Generate Neumann regions used for an example
%   Returns a segments object composed of multiple star-shaped regions
%   centered at specific complex coordinates, using the mystar constructor
%   from the 'segments' class (as in mpspack).

% Define centers of the regions (complex numbers)
cs = [
    0.720 + 0.353i;
    0.320 + 0.420i;
    0.540 + 0.508i;
    0.749 + 0.704i;
    0.408 + 0.725i;
    0.130 + 0.276i;
    0.133 + 0.907i;
    0.369 + 0.169i
];

% Define radii and oscillation numbers for each region
rs = [
    0.126;
    0.081;
    0.082;
    0.135;
    0.118;
    0.108;
    0.071;
    0.071
];

oscil = [6; 5; 3; 3; 6; 5; 6; 7];

% Number of regions
M = numel(cs);

% Create cell array of segment objects
segs = cell(M, 1);
for i = 1:M
    segs{i} = segments.mystar(cs(i), rs(i), oscil(i));
end

% Construct the segments object and set properties
ss = segments(segs);
ss.centers(cs);
ss.radiuss(rs);
ss.oscil(oscil);

end

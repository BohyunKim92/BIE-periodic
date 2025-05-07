function [ss] = circles(cs, rs)
%CIRCLES Generate a segments object composed of circular segments.
%   ss = CIRCLES(cs, rs) creates a segments object consisting of circles
%   centered at complex points in cs with corresponding radii in rs.
%
%   Inputs:
%       cs - complex array of centers for the circles
%       rs - real array of radii for the circles
%
%   Output:
%       ss - segments object containing all circular components

M = numel(cs);
segs = cell(M, 1);
oscil = ones(M, 1);  % oscillation parameter = 1 for circles

for i = 1:M
    segs{i} = segments.mycircle(cs(i), rs(i));
end

ss = segments(segs);
ss.centers(cs);
ss.radiuss(rs);
ss.oscil(oscil);

end

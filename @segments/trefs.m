function [ss] = trefs(cs, rs)
%TREFS Generate segments object of trefoil-shaped boundaries
%   [ss] = TREFS(cs, rs) returns a 'segments' object composed of
%   trefoil-shaped domains, where `cs` are complex centers and `rs` are radii.
%
%   Each trefoil has a fixed oscillation number of 3.

% Number of trefoil shapes
M = numel(cs);

% Initialize cell array for segment objects and set fixed oscillation
segs = cell(M, 1);
oscil = 3 * ones(M, 1);

% Generate each trefoil segment
for i = 1:M
    segs{i} = segments.mytref(cs(i), rs(i));
end

% Create segments object and assign attributes
ss = segments(segs);
ss.centers(cs);
ss.radiuss(rs);
ss.oscil(oscil);

end

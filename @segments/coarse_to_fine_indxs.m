function indxs = coarse_to_fine_indxs(finess, coarsess)
%COARSE_TO_FINE_INDXS Maps coarse discretization indices to fine ones.
%   indxs = COARSE_TO_FINE_INDXS(finess, coarsess) returns a mapping from 
%   the coarse to fine indices, assuming that the refinement factor between 
%   fine and coarse is a power of 5: finess.Ns ./ coarsess.Ns = 5^p.
%
%   Inputs:
%       finess   - Fine discretization structure with fields:
%                    - Ns: number of nodes per segment
%                    - indxs: index arrays for fine segments
%       coarsess - Coarse discretization structure with fields:
%                    - Ns: number of nodes per segment
%                    - indxs: index arrays for coarse segments
%                    - tN: total number of nodes
%
%   Output:
%       indxs - Vector of fine-level indices corresponding to coarse nodes

M = numel(finess.Ns);
factors = finess.Ns ./ coarsess.Ns;
pows = round(log10(factors) ./ log10(5));  % Determine the power p such that 5^p = factor

indxs = zeros(coarsess.tN, 1);

for i = 1:M
    stride = 5^pows(i);
    start_idx = ceil(stride / 2);  % Centered sample
    temp = start_idx : stride : finess.Ns(i);
    indxs(coarsess.indxs{i}) = finess.indxs{i}(temp);
end

end

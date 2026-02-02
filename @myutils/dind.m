function [inddiag] = dind(A)
% Given a square matrix A, get the index of diagonal of the matrix.
N = size(A,1); % size of the square matrix
inddiag = sub2ind(size(A),1:N,1:N);
end


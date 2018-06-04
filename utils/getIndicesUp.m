function Iu = getIndicesUp(N)
% Compute the linear indices corresponding to the upper triangular part of 
% a square matrix of size N.
%-------------------------------------------------------------------------%
% Input:
% > N : matrix size
%
% Output:
% < Iu : indices corresponding to the upper triangular part of the matrix
%
%-------------------------------------------------------------------------%
%% 
% [15/03/2018], P.-A. Thouvenin.
%-------------------------------------------------------------------------%

Iu = zeros(N*(N-1)/2,1);
k = 1;
for j = 2:N
    Iu(k:k+j-2) = (j-1)*N + (1:j-1)';
    k = k + j - 1;
end

end


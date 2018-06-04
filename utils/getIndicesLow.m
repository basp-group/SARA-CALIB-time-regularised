function Il = getIndicesLow(N)
% Compute the linear indices corresponding to the lower triangular part of 
% a square matrix of size N.
%-------------------------------------------------------------------------%
% Input:
% > N : matrix size
%
% Output:
% < Il : indices corresponding to the lower triangular part of the matrix
%
%-------------------------------------------------------------------------%
%% 
% [15/03/2018], P.-A. Thouvenin.
%-------------------------------------------------------------------------%

Il = zeros(N*(N-1)/2,1);
k = 1;
for j = 1:N-1
    Il(k:k+N-j-1) = (j-1)*N + (j+1:N)';
    k = k + N - j;
end

end


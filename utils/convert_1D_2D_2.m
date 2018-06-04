function c = convert_1D_2D_2(n,K)
% Converts 1D linear indices to 2D indices (matrix indexing).
%-------------------------------------------------------------------------%
% Input:
% > n : vector containing the linear indices of interest [N, 1]
% > K : matrix size [1, 2]
%
% Note: n(i) \in {1, ..., K(1)*K(2)}
%
% Output:
% < n : matrix containing the corresponding 2D indices [N, 2]
% 
% Note: c(i, :) \in \{-K(1)/2, ..., K(1)/2-1}x{-K(2)/2, ..., K(2)/2-1}
%
%-------------------------------------------------------------------------%
%% 
% [15/03/2018], P.-A. Thouvenin.
%-------------------------------------------------------------------------%

d = floor(n/K(1)) ;
c = [n - (d*K(1)) - floor(K(1)/2) - 1, d - floor(K(2)/2)] ;

end
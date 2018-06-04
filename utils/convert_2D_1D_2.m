function n = convert_2D_1D_2(c, K)
% Converts 2D indices (matrix indexing) to linear indices.
%-------------------------------------------------------------------------%
% Input:
% > c : matrix containing the 2D indices of interest [N, 2]
% > K : matrix size [1, 2]
%
% Note: c(i, :) \in \{-K(1)/2, ..., K(1)/2-1}x{-K(2)/2, ..., K(2)/2-1}
%
% Output:
% < n : corresponding linear indices [N, 1]
% 
% Note: n(i) \in {1, ..., K(1)*K(2)}
%
%-------------------------------------------------------------------------%
%% 
% [15/03/2018], P.-A. Thouvenin.
%-------------------------------------------------------------------------%

n = (c(:,2) + floor(K(2)/2))*K(1) + c(:,1) + floor(K(1)/2) +1;

end
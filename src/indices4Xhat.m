function id = indices4Xhat(S, J, K, om_t)
% Compute the indices of the elements to be extracted from xhat to .
%-------------------------------------------------------------------------%
% Input:
% > X     : part of xhat [S2, S2, J2, na - 1]
% < V     : part of the spatial NUFFT gridding kernels with redundancy (similar to Y) [J2, na - 1]
% > D1_at : D1(:, :, t), where the column a has been removed 
%           D1 : [S2, na, T] -> D1_at : [S2, na-1]
% > J     : size of the gridding kernels
% > S     : size of the DDE support (spatial Fourier domain)
% > id_nnz: position of missing measurements (temporal dimension)
%
% Output:
% < H1a   : linear operator Ha [S2, na - 1]
%
%-------------------------------------------------------------------------%
%% 
% [15/03/2018], P.-A. Thouvenin.
%-------------------------------------------------------------------------%
%%

% om_t: [na-1, 2]

% Auxiliary variables
S2 = S^2; % assuming square DDE supports
J2 = J^2;
n = size(om_t, 1); % n = na - 1

%% Generate the indices to select the columns (area centered in the 0 frequency) (total area [S,S])
ci = nufft_offset(0, S, K(1)); % [1, 1] position of the leftmost element in the J nbr
cj = nufft_offset(0, S, K(2)); % [1, 1] position of the leftmost element in the J nbr
kdy = mod(bsxfun(@plus, (1:S)', ci.'), K(1)) + 1; % [S 1]
kdx = mod(bsxfun(@plus, (1:S)', cj.'), K(2)) + 1; % [S 1] 
kdx = (kdx - 1)*K(1); % -> pre-convert these indices into offsets
ll0 = reshape(bsxfun(@plus, reshape(kdy, [S, 1]), reshape(kdx, [1, S])), [S2, 1]); % [S2 1]
ll0 = fftshiftId2D(K(1),K(2),ll0(:)); % [S2, 1]

%% Compute the indices for the elements to be fetched in x_hat
ciJ = nufft_offset(om_t(:, 1), J, K(1)); % [n, 1] position of the center
cjJ = nufft_offset(om_t(:, 2), J, K(2)); % [n, 1] position of the center
kdyJ = mod(bsxfun(@plus, (1:J)', ciJ.'), K(1)) + 1; % [J n]
kdxJ = mod(bsxfun(@plus, (1:J)', cjJ.'), K(2)) + 1; % [J n] 
kdxJ = (kdxJ - 1)*K(1); % -> pre-convert these indices into offsets
llJ = reshape(bsxfun(@plus, reshape(kdyJ, [J, 1, n]), reshape(kdxJ, [1, J, n])), [J2, n]); % [J2 n]
llJ = fftshiftId2D(K(1), K(2),llJ(:)); % [J2*n, 1]
dJ = convert_1D_2D_2(llJ, K); % [J2*n, 2]

d = convert_1D_2D_2(ll0, K); % [S2, 2]
d = reshape(bsxfun(@plus, reshape(d, [S2, 1, 2]), reshape(d, [1, S2, 2])), [S2^2, 1, 2]); % [S2,S2,2]
d = bsxfun(@plus, d, reshape(dJ, [1, J2*n, 2])); % [S2^2, J2*n, 2]
id = convert_2D_1D_2(reshape(d, [(S2^2)*J2*n, 2]), K); % [(S2^2)*J2*n, 1]

end

function [G, ll, v, V] = createGnufft_T2(D1, D2, om, K, S, J, W)
% Create the gridding matrix G including the DDE kernels.
%-------------------------------------------------------------------------%
% Input:
% > D1, D2 : DDEs kernels for a single time instant [S2, na]
% > om     : normalized u-v frequencies [M, 2] 
% > K      : size of the spatial Fourier space [1, 2] 
% > S      : size of the DDE kernels (in the spatial Fourier domain)
% > J      : size of the gridding kernels
% > W      : values of the gridding kernels [J2, M]
%
% Output:
% < G  : gridding matrix [M, prod(K)] (sparse matrix)
% < ll : position of the nonzero values of G
% < v  : convolutions between the DDE kernels [Q2, M]
% < V  : values contained in G [Q2, J2, M]
%
%-------------------------------------------------------------------------%
%% 
% [08/12/2017], P.-A. Thouvenin.
%-------------------------------------------------------------------------%
%%
Q  = 2*S - 1;             % size of the kernels after convolution
Q2 = Q^2;
J2 = J^2;

Qprime = floor(Q/2);      % Q is always odd (one possibility only)
tmp1 = (-Qprime:Qprime).';
tmp1 = tmp1(:, ones(Q,1));

[~, na] = size(D1);    % [S2, na] number of antennas at time t
M_true  = na*(na-1)/2; % number of acquisitions -> check value...
M = size(om, 1);       % M = M_true if all the measurements are present, M < M_true if some of the data have been flagged
v = zeros(Q^2,M_true); % convolution values (stored in column for each pair)
q = 0;                 % global counter

%% Perform 2D convolutions and gridding using D1 and D2 (to be possibly performed in parallel)
for alpha = 1:na-1
    for beta = alpha+1:na % modify the double loop to exclusively select the appropriate elements, apply nonzeros on W
        % 2D convolutions
        q = q+1;
        v(:,q) = reshape(conv2(rot90(reshape(D1(:,alpha),[S,S]),2),reshape(conj(D2(:,beta)),[S,S])), [Q^2,1]); % only select the appropriate entries...
    end
end

% Generate indices in the sparse G matrix
if rem(J,2) > 0 % odd
   c0 = round(om.*K/(2*pi)) - (J+1)/2; % [M, 2]
else
   c0 = floor(om.*K/(2*pi)) - J/2; 
end
kdy = bsxfun(@plus, (1:J).', c0(:,1).'); % [J M]
kdx = bsxfun(@plus, (1:J).', c0(:,2).'); % [J M]
ii = mod(bsxfun(@plus, tmp1(:), reshape(kdy, [1,J,M])), K(1)) + 1; % [Q2, J, n] % row indices of the elements within each area, 
                                                                   % whose leftmost element row indices are given above
jj = mod(bsxfun(@plus, reshape(tmp1.', [Q2,1]), reshape(kdx, [1,J,M])), K(2)) + 1; % [Q2, J, M] % column indices ...
ll = reshape(bsxfun(@plus, reshape((jj-1)*K(1), [Q2,1,J,M]), reshape(ii, [Q2,J,1,M])), [Q2,J2,M]);

% Duplicate values to have all the convolutions centered in the different elements
V = bsxfun(@times, reshape(v, [Q2, 1, M_true]), reshape(W, [1, J2, M_true])); % [Q2, J2, M]   
V(isnan(V(:))) = []; % [J2*M, 1] there are zeros in W at the positions where the 
                 % measurements are missing, right size once the zeros are 
                 % filtered

% Generate row indices (within G) [remark: Jt^2 nz elements per row]
kk = bsxfun(@times, ones(J2*Q2,1), 1:M); % [J2*Q2, M]
% try
G = sparse(kk(:), ll(:), V(:), M, K(1)*K(2));  % [M, prod(K)]
% catch 
%     keyboard
% end
end

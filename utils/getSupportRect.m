function [I,Ic] = getSupportRect(k, S, N1, N2)
% Generate the indices corresponding to the pixels contained in the 
% frequency region defined by the support S, centered in the spatial
% frequencies given in K.
%
% Input
% > k     : spatial frequencies, in {1,...,N1N2}, [Na, 1]
% > S     : support size (assume square support)
% > N1,N2 : image dimensions.
%
% Output
% < I  : linear indices of the pixel contained in the frequency region of 
%        size S, centered in the elements of K (linear indices) [S^2*Na, 1]
% < Ic : indices corresponding to the entries on I [S^2*Na, 2].
%-------------------------------------------------------------------------%
%%
% Remark:
% entries to be filtered out with ~isnan(), output format to be modified;
%-------------------------------------------------------------------------%
%%
% Code : Pierre-Antoine Thouvenin, October 9th 2017.
% Correction : [14/11/2017]
%%
% Auxiliary variables (to handle the cases where S is even or odd)
r = rem(S,2) - 1;
S_prime = floor(S/2);
Na = size(k,1); % number of frequencies to be processed simultaneously

% Convert the linear indices of the centers to 2D indices (ind2sub)
[K1,K2] = ind2sub([N1,N2],k);
K = [K1,K2];

% Center the elements in [1,N]^2
% K_prime = K + floor(N/2) + 1; % assuming the frequencies are already
% centered as appropriate

% Generate all the couples (3D)
tmp = (-S_prime:(S_prime + r))'; % row indices defining the elements in S
tmp2 = tmp(:,ones(S,1));         % column idices defining the elements in S
% keyboard
n_c = [reshape(tmp2, [S^2,1]), reshape(tmp2', [S^2,1])]; 
Ic = bsxfun(@plus, n_c, reshape(K', [1, 2, size(K, 1)])); % [S^2,2,Na]
Ic = reshape(permute(Ic, [1,3,2]), [S^2*Na,2]);           % [S^2*Na,2]

% Filter out invalid positions (i.e., out of the scope of the image)
id1 = bsxfun(@or, Ic(:,1) < 1, Ic(:,1) > N1);
id2 = bsxfun(@or, Ic(:,2) < 1, Ic(:,2) > N2);
Ic(id1,:) = NaN; % [S^2*Na, 2]
Ic(id2,:) = NaN; % [S^2*Na, 2]

% Convert to linear indices
I = (Ic(:,2) - 1)*N1 + Ic(:,1);

end

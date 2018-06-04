function D = computeD(U, F, Gt, scale, T)
% Compute the DDE kernels in the temporal domain.
%-------------------------------------------------------------------------%
% Input:
% > U     : DDEs in the Fourier domain (spatial and temporal) [S2, na, P]
% > F     : size of the temporal Fourier domain
% > Gt    : temporal gridding matrix (unused, for possible extension to temporal NUFFT)
% > scale : scaling factor (unused, for possible extension to temporal NUFFT)
% > T     : number of time instants
%
% Output:
% < D     : DDEs in the temporal domain [S2, na, T]
%
%-------------------------------------------------------------------------%
%% 
% [21/01/2018] Debug and check values, possibly add 1D NUFFT (temporal dimension)
% [15/03/2018], P.-A. Thouvenin.
%-------------------------------------------------------------------------%
%%
[S2, n, P] = size(U);
% T = size(Gt, 1);
% Compute D1 (or D2) for the 2D convolutions to be performed later
% D = Gt*so_ifft(reshape(permute(U, [3, 1, 2]), [P, S2*n]), T, F, scale); % [T, n*S2], transform (fft or ifft) performed along the dimension 1 for U1/U2
% % keyboard
% D = so_ifft(reshape(permute(U, [3, 1, 2]), [P, S2*n]), T, F, scale);
% D = D(1:T, :);
% D = reshape(D.', [S2, n, T]);

% Third version
x = zeros(T, S2*n);
c = floor(T/2) + 1;
p = floor(P/2);
x(c-p:c+p, :) = reshape(permute(U, [3, 1, 2]), [P, S2*n]); % [T, S2*n]
D = reshape(T*ifft(ifftshift(x, 1), T, 1).', [S2, n, T]); % sqrt(T)

end

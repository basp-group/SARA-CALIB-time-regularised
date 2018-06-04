function Ua = computeUa(Da, F, T, P, Gt, scale)
% Compute the DDE kernel in the temporal Fourier domain for a single antenna.
%-------------------------------------------------------------------------%
% Input:
% < Da     : DDEs in the temporal domain [S2, T]
% > F     : size of the temporal Fourier domain
% > T     : number of time instants
% > P     : size of the DDE support in the temporal Fourier domain
% > Gt    : temporal gridding matrix (unused, for possible extension to temporal NUFFT)
% > scale : scaling factor (unused, for possible extension to temporal NUFFT)
% > id_a  : position of possibly missing measurements in the temporal
%           domain
%
% Output:
% < Ua    : DDEs in the Fourier domain for antenna a (spatial and temporal) [S2, P]
%
%-------------------------------------------------------------------------%
%% 
% [21/01/2018] Debug and check values, possibly add 1D NUFFT (temporal dimension)
% [15/03/2018], P.-A. Thouvenin.
%-------------------------------------------------------------------------%
%%
% [S2, T] = size(Da);
% [S2, P] = size(Ua);
% Compute D1 (or D2) for the 2D convolutions to be performed later
% Ua = so_ifft_adj((Da*conj(Gt)).', P, F, scale).'; % [S2, P] % zero entries are filtered from the transform 
%%- Ua = so_ifft_adj((Da(:, id_t)*conj(Gt(id_t,:))).', P, F, scale).';

% Ua = so_ifft_adj(Da.', T, P, scale).'; % [S2, P]

% Corrected version
% c = floor(T/2) + 1;
% p = floor(P/2);
% y = fftshift(fft(Da.', T, 1), 1); % [T, S2]
% Ua = y(c-p:c+p, :).'; % [S2, P]

% Corrected version
c = floor(F/2) + 1;
p = floor(P/2);
y = fftshift(fft(Da.', F, 1), 1); % [T, S2] /sqrt(F)
Ua = y(c-p:c+p, :).'; % [S2, P]

end

function Da = computeDa(Ua, F, Gt, scale, T, id_a)
% Compute the DDE kernel in the temporal domain for a single antenna.
%-------------------------------------------------------------------------%
% Input:
% > Ua    : DDEs in the Fourier domain for antenna a (spatial and temporal) [S2, P]
% > F     : size of the temporal Fourier domain
% > Gt    : temporal gridding matrix (unused, for possible extension to temporal NUFFT)
% > scale : scaling factor (unused, for possible extension to temporal NUFFT)
% > T     : number of time instants
% > id_a  : position of possibly missing measurements in the temporal
%           domain
%
% Output:
% < Da     : DDEs in the temporal domain [S2, T]
%
%-------------------------------------------------------------------------%
%% 
% [21/01/2018] Debug and check values, possibly add 1D NUFFT (temporal dimension)
% [15/03/2018], P.-A. Thouvenin.
%-------------------------------------------------------------------------%
%%
% Version with Temporal NUIFFT
% D: [T, na*S2]
% U: [P, na*S2]
% [S2, P] = size(Ua);
% Compute D1 (or D2) for the 2D convolutions to be performed later
% Da = (Gt(id_t,:)*so_ifft(Ua.', F, scale)).'; % [S2, T'] % zero entries are filtered from the transform
% Da = (Gt*so_ifft(Ua.', F, scale)).'; % [S2, T]
% Da = so_ifft(Ua.', T, F, scale);
% Da = Da(1:T, :).'; % [S2, T]

% % Corrected version
% [S2, P] = size(Ua);
% x = zeros(T, S2);
% c = floor(T/2) + 1;
% p = floor(P/2);
% x(c-p:c+p, :) = Ua.';
% Da = T*ifft(ifftshift(x, 1), T, 1).'; % [S2, T]

% Test with additional 0-padding
[S2, P] = size(Ua);
x = zeros(F, S2);
c = floor(F/2) + 1;
p = floor(P/2);
x(c-p:c+p, :) = Ua.';
Da = sqrt(F)*ifft(ifftshift(x, 1), F, 1).'; % [S2, T] % F*
Da = Da(:, 1:T);
Da(id_a) = 0; % set to 0 the time instants corresponding to missing measurements

end

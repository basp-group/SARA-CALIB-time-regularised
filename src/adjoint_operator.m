function x = adjoint_operator(y, Ha, F, T, Gt, scale, P)
% Compute the adjoint of the linear operator involved in the update of the
% DDEs for a single antenna (na antennas in total).
%-------------------------------------------------------------------------%
% Input:
% > y     : input vector [na - 1, T]
% > Ha    : linear operator Ha [na-1, S2, T]
% > F     : size of the temporal Fourier domain
% > T     : number of time instants
% > Gt    : temporal gridding matrix (unused, for possible extension to temporal NUFFT)
% > scale : scaling factor (unused, for possible extension to temporal NUFFT)
% > P     : size of the DDE support (temporal Fourier domain)
%
% Output:
% < x     : output vector [S2, P]
%
%-------------------------------------------------------------------------%
%% 
% [15/03/2018], P.-A. Thouvenin.
%-------------------------------------------------------------------------%
%%
[n, S2, T] = size(Ha); % n = na - 1
x = computeUa(reshape(sum(bsxfun(@times, permute(conj(Ha), [2, 1, 3]), ...
    reshape(y, [1, n, T])), 2), [S2, T]), F, T, P, Gt, scale); % [S2, P]

end

function y = direct_operator(x, Ha, F, T, Gt, scale, id_a)
% Compute the linear operator involved in the update of the
% DDEs for a single antenna (na antennas in total).
%-------------------------------------------------------------------------%
% Input:
% > x     : input vector [S2, P]
% > Ha    : linear operator Ha [na-1, S2, T]
% > F     : size of the temporal Fourier domain
% > T     : number of time instants
% > Gt    : temporal gridding matrix (unused, for possible extension to temporal NUFFT)
% > scale : scaling factor (unused, for possible extension to temporal NUFFT)
% > id_a  : vector containing the indices of possibly missing measurements
%           (in the temporal domain)
%   P     : size of the DDE support (temporal Fourier domain)
%
% Output:
% < x     : output vector [na - 1, T] 
%
%-------------------------------------------------------------------------%
%% 
% [15/03/2018], P.-A. Thouvenin.
%-------------------------------------------------------------------------%

[~, S2, T] = size(Ha); % T2 <= T (missing entires removed from the total number of time instants T)
y = squeeze(sum(bsxfun(@times, Ha, reshape(computeDa(x, F, Gt, scale, T, id_a), [1, S2, T])), 2)); % [na-1, T2]
% leave 0 in H for missing entries
end

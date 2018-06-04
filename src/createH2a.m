function H2a = createH2a(X, V, D1_at, J, S, id_nnz)
% Compute the linear operator Ha involved in the update of the second
% DDE term for a single antenna a (na antennas in total).
%-------------------------------------------------------------------------%
% Input:
% > X     : part of xhat [S2*S2*J2*(na - 1), 1]
% > V     : part of the spatial NUFFT gridding kernels with redundancy (similar to Y) [J2, na - 1]
% > D1_at : D1(:, :, t), where the column a has been removed 
%           D1 : [S2, na, T] -> D1_at : [S2, na-1]
% > J     : size of the gridding kernels
% > S     : size of the DDE support (spatial Fourier domain)
% > id_nnz: position of missing measurements (temporal dimension)
%
% Output:
% < H2a   : linear operator Ha [S2, na - 1]
%
%-------------------------------------------------------------------------%
%% 
% [15/03/2018], P.-A. Thouvenin.
%-------------------------------------------------------------------------%
%%

S2 = S^2;
J2 = J^2;
n = length(id_nnz); % n = na - 1 if all the measurements are present, n < na - 1 otherwise. 

H2a = zeros(n, S2);
for q = 1:n
    Xq = sum(bsxfun(@times, reshape(V(:,id_nnz(q)), [1, 1, J2]), reshape(X((q-1)*S2^2*J2+1:q*S2^2*J2), [S2, S2, J2])), 3); % to be verified
    H2a(q,:) = flipud(D1_at(:,q)).'*Xq; % select the appropriate portion of D2_at to account for missing measurements
end

end


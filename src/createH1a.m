function H1a = createH1a(X, V, D2_at, J, S, id_nnz)
% Compute the linear operator Ha involved in the update of the first
% DDE term for a single antenna a (na antennas in total).
%-------------------------------------------------------------------------%
% Input:
% > X     : part of xhat [S2*S2*J2*(na - 1), 1]
% > V     : part of the spatial NUFFT gridding kernels with redundancy (similar to Y) [J2, na - 1]
% > D2_at : D2(:, :, t), where the column a has been removed 
%           D2 : [S2, na, T] -> D2_at : [S2, na-1]
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
S2 = S^2;
J2 = J^2;
n = length(id_nnz); % n = na - 1 if all the measurements are present, n < na - 1 otherwise. 

% H1a = zeros(S2, n);
% for q = 1:n
%     Xq = sum(bsxfun(@times, reshape(V(:,id_nnz(q)), [1, 1, J2]), reshape(X((q-1)*S2^2*J2+1:q*S2^2*J2), [S2, S2, J2])), 3); % to be verified
% %     Xq = sum(bsxfun(@times, reshape(V(:,q), [1, 1, J2]), reshape(X{q}, [S2, S2, J2])), 3); % to be verified, currently error with reshape
%     H1a(:,q) = Xq*conj(D2_at(:,q)); % select the appropriate portion of D2_at to account for missing measurements
%     % H1a(:,q) = Xq*conj(D2_at(:,id_nnz(q)));
% end

% V: [J2, n]
Xq = reshape(sum(bsxfun(@times, reshape(V, [J2, 1, 1, n]), permute(reshape(X, [S2, S2, J2, n]), [3, 1, 2, 4])), 1), [S2, S2, n]); % [S2, S2, n]
H1a = reshape(sum(bsxfun(@times, permute(Xq, [2, 1, 3]), reshape(conj(D2_at), [S2, 1, n])), 1), [S2, n]); % [S2, n]

end

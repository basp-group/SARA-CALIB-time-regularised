function Ua2 = updateU2a(Ya, Ua1, Ua2, H2a, Gt, scale_t, F, T, JU2o, nuo, mu, lambda_scale, theta_minoR, theta_maxoR, theta_minoI, theta_maxoI)
% Update the DDE term 2 related to antenna a (PALM step).
%-------------------------------------------------------------------------%
% Input:
% > Ya           : column a of the data matrix Y with redundancy [na-1, T]
% > Ua2          : DDEs related to antenna a [S2, P]
% > Ua1          : DDEs related to antenna a [S2, P]
% > H2a          : linear operator involved in the update step [S2, na - 1]
% > Gt           : (unused) gridding matrix for temporal NUFFT
% > scale_t      : (unused) scaling coeff. for temporal NUFFT
% > F            : size of the temporal Fourier space
% > T            : number of time instants
% > JU2o         : maximum number of iterations
% > nuo          : regularization parameter (distance between U1a and U2a)
% > mu           : regularization parameter (distance between U1a and Phi)
% > lambda_scale : = 1/prod(K), K size of the spatial Fourier space
% > theta_minoR, theta_maxoR : box constraints on the real part of U1a
% > theta_minoI, theta_maxoI : box constraints on the imaginary part of U1a
%
% Output:
% < Ua2          : DDEs related to antenna a [S2, P]
%
%-------------------------------------------------------------------------%
%% 
% [15/03/2018], P.-A. Thouvenin.
%-------------------------------------------------------------------------%
%%

% U1: [S2, na, P] -> D1: [S2, na, T] 
% U2: [na, S2, P] -> D2: [na, S2, T] 
% Y: [na, na, T] -> Ya: [na-1, T]
% Ua1: [S2, P]

[S2, P] = size(Ua2);
u1a = conj(Ua1); % [S2, P]
u2a = conj(Ua2); % [S2, P]
Phi = sparse(floor(S2/2)+1, floor(P/2)+1, 1, S2, P);
id_a = find(Ya == 0); % indices of the missing measurements (represented by 0) 

% step size                    
A = @(x) direct_operator(fliplr(x), H2a, F, T, Gt, scale_t, id_a); % do not forget the fliplr (property F(x*)[n] = {F(x)*}[-n]) 
At = @(x) fliplr(adjoint_operator(x, H2a, F, T, Gt, scale_t, P));
Lips_temp = mu + nuo + 2*lambda_scale*pow_method(@(x) A(x), @(x) At(x), size(u2a))/F;
gamma1 = 1.9/Lips_temp;
    
% ----------------------------------------------------
% Iterations
% ----------------------------------------------------
for q = 1:JU2o
    % gradient step
    grad = 2*lambda_scale*At(squeeze(A(u2a)) - Ya)/F + nuo*(u2a - u1a) + mu*(u2a - Phi); % [S2, P], 0 are preserved using the structure of H2_a
    g = u2a - gamma1*grad;
    % proximity step
    vr = min(max(real(g), theta_minoR), theta_maxoR);
    vi = min(max(imag(g), theta_minoI), theta_maxoI);
    u2a = vr + 1i*vi;
end
% ----------------------------------------------------
Ua2 = conj(u2a); % [S2, P]

end


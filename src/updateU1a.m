function Ua1 = updateU1a(Ya, Ua1, Ua2, H1a, Gt, scale_t, F, T, JU1o, nuo, mu, lambda_scale, theta_minoR, theta_maxoR, theta_minoI, theta_maxoI)
% Update the DDE term 1 related to antenna a (PALM step).
%-------------------------------------------------------------------------%
% Input:
% > Ya           : column a of the data matrix Y with redundancy [na-1, T]
% > Ua1          : DDEs related to antenna a [S2, P]
% > Ua2          : DDEs related to antenna a [S2, P]
% > H1a          : linear operator involved in the update step [S2, na - 1]
% > Gt           : (unused) gridding matrix for temporal NUFFT
% > scale_t      : (unused) scaling coeff. for temporal NUFFT
% > F            : size of the temporal Fourier space
% > T            : number of time instants
% > JU1o         : maximum number of iterations
% > nuo          : regularization parameter (distance between U1a and U2a)
% > mu           : regularization parameter (distance between U1a and Phi)
% > lambda_scale : = 1/prod(K), K size of the spatial Fourier space
% > theta_minoR, theta_maxoR : box constraints on the real part of U1a
% > theta_minoI, theta_maxoI : box constraints on the imaginary part of U1a
%
% Output:
% < Ua1          : DDEs related to antenna a [S2, P]
%
%-------------------------------------------------------------------------%
%% 
% [15/03/2018], P.-A. Thouvenin.
%-------------------------------------------------------------------------%
%%

% U1: [S2, na, P] -> D1: [S2, na, T] 
% U2: [S2, na, P] -> D2: [S2, na, T] 
% Y: [na, na, T] -> Ya: [na-1, T], with 0 for missing measurements
% Ua1: [S2, P]
% H1a: [na-1, S2, T]

[S2, P] = size(Ua1);
u1a = flipud(Ua1); % [S2, P]
u2a = flipud(Ua2); % [S2, P]
Phi = sparse(floor(S2/2)+1, floor(P/2)+1, 1, S2, P);
id_a = find(Ya == 0); % indices of the missing measurements (represented by 0) 

% step size                              
A = @(x) direct_operator(x, H1a, F, T, Gt, scale_t, id_a);
At = @(x) squeeze(adjoint_operator(x, H1a, F, T, Gt, scale_t, P));
Lips_temp = mu + nuo + 2*lambda_scale*pow_method(@(x) A(x), @(x) At(x), size(u1a)); % /F
gamma1 = 1.9/Lips_temp;

% ----------------------------------------------------
% Iterations
% ----------------------------------------------------
for q = 1:JU1o
    % gradient step
    grad = 2*lambda_scale*At(squeeze(A(u1a)) - Ya) + nuo*(u1a - u2a) + mu*(u1a - Phi); % /F [S2, P], 0 are preserved using the structure of H1_a
    g = u1a - gamma1*grad;
    % proximity step
    vr = min(max(real(g), theta_minoR), theta_maxoR);
    vi = min(max(imag(g), theta_minoI), theta_maxoI);
    u1a = vr + 1i*vi;
end
% ----------------------------------------------------
Ua1 = flipud(u1a); % [S2, P]

end


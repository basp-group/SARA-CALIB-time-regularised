function Ut1 = updateU1t(Yt, Ut1, Ut2, Ht1, na_t, JU1o, nuo, mu, lambda_scale, theta_minoR, theta_maxoR, theta_minoI, theta_maxoI)
% Update the DDE term 1 for a given time instant t (PALM step).
%-------------------------------------------------------------------------%
% Input:
% > Yt           : data matrix with redundancy at time instant t [na, na]
% > Ut1          : DDEs at time instant t [S2, na]
% > Ut2          : DDEs at time instant t [S2, na]
% > Ht1          : linear operator involved in the update [na-1, S2, na]
% > na_t         : list of valid antennas (not masked)
% > JU1o         : maximum number of iterations
% > nuo          : regularization parameter (distance between U1a and U2a)
% > mu           : regularization parameter (distance between U1a and Phi)
% > lambda_scale : = 1/prod(K), K size of the spatial Fourier space
% > theta_minoR, theta_maxoR : box constraints on the real part of U1a
% > theta_minoI, theta_maxoI : box constraints on the imaginary part of U1a
%
% Output:
% < Ut1          : DDEs related to antenna a [S2, na]
%
%-------------------------------------------------------------------------%
%% 
% [15/03/2018], P.-A. Thouvenin.
%-------------------------------------------------------------------------%
%%

% U2t: [S2, na]
% U1t: [S2, na]
% Yt: [na, na]
% H1t: [na-1, S2, na]
id = true(na_t,1);
S2 = size(Ut1, 1);
phi = sparse(floor(S2/2)+1, 1, 1, S2, 1);

for a = 1:na_t
    id(a) = false;
    
    % construction of data vector
    Yt_alpha = nonzeros(Yt(a,id));  % transpose of the rows of Y /!\ nonzeros systematically returns a column vector
    u1t = flipud(Ut1(:,a));
    u2t = flipud(Ut2(:,a));
    
    % step size
    Lips_temp = mu + nuo + 2*lambda_scale*pow_method(@(x) Ht1(:, :, a)*x, @(x) Ht1(:, :, a)'*x, size(u1t));
    gamma1 = 1.9/Lips_temp ;

    % ----------------------------------------------------
    % Iterations
    % ----------------------------------------------------
    for q = 1:JU1o
        % gradient step
        grad1 = 2*lambda_scale*Ht1(:, :, a)'*(Ht1(:, :, a)*u1t - Yt_alpha);
        grad2 = nuo*(u1t - u2t) + mu*(u1t - phi);
        grad = grad1 + grad2;
        g = u1t - gamma1*grad;
        % proximity step
        vr = min(max(real(g), theta_minoR), theta_maxoR);
        vi = min(max(imag(g), theta_minoI), theta_maxoI);
        u1t = vr + 1i*vi;
    end
    % ----------------------------------------------------
    Ut1(:,a) = flipud(u1t);
    id(a) = true;
end

end


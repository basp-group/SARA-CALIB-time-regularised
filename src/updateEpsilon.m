function epsilon = updateEpsilon(y, x0, epsilon, Jeps, tol_x, nit_tot, proxEpsilon, B_tmp, Bt_tmp, lambda_scale)
% Update the image parameter epsilon (PALM descent steps).
%-------------------------------------------------------------------------%
% Input:
% > y            : data vector [M, 1]
% > x0           : reference image term [N]
% > epsilon      : initial epsilon term [N]
% > Jeps         : maximum number of iterations
% > tol_x        : stopping criterion
% > proxEpsilon  : lambda function for the proximal operator
% > B_tmp        : lambda function representing the direct measurement operator
% > Bt_tmp       : adjoint of B_tmp (lambda function)
% > lambda_scale : rescaling term (used everywhere in the problem)
%                  lambda_scale = 1/prod(K), K dimension of the spatial
%                  Fourier space
%
% Output:
% < epsilon   : updated epsilon term
%
%-------------------------------------------------------------------------%
%% 
% [15/03/2018], P.-A. Thouvenin.
%-------------------------------------------------------------------------%
%%

Ax = 2*lambda_scale*pow_method(B_tmp, Bt_tmp, size(epsilon));
ytmp = y - B_tmp(x0);

for q = 1:Jeps
    eps_old = epsilon;
    grad_eps = 2*lambda_scale*Bt_tmp(B_tmp(epsilon) - ytmp);
    eps_tmp = epsilon - (1.9/Ax)*real(grad_eps); % to be verified
    epsilon = proxEpsilon(eps_tmp, Ax);
     
    % stopping criterion (imaging step)
    if (q>10) && (norm(epsilon - eps_old, 'fro')/norm(epsilon, 'fro') < tol_x)
        disp(['x: stopping criterion reached, glob. iteration = ', num2str(nit_tot)])
        disp(['x: stopping criterion reached, inner iteration = ', num2str(q)])
        break
    end
end

end

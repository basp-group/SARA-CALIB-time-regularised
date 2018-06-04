function [sol, norm_l1] = solver_prox_L1_full_image(x, lambda, param)
% PROJ_L1 - Proximal operator with L1 norm
%
% sol = solver_prox_L1_full_image(x, lambda, param) solves:
%
%   min_{z} 0.5*||x - z||_2^2 + lambda * ||Psit (xA + x)||_1
%
% The input argument param contains the following fields:
%
%   - Psit: Sparsifying transform (default: Id).
%
%   - Psi: Adjoint of Psit (default: Id).
%
%   - tight: 1 if Psit is a tight frame or 0 if not (default = 1)
%
%   - nu: bound on the norm of the operator A, i.e.
%       ||Psi x||^2 <= nu * ||x||^2 (default: 1)
%
%   - max_iter: max. nb. of iterations (default: 200).
%
%   - rel_obj: minimum relative change of the objective value (default:
                                                               %   1e-4)
%       The algorithm stops if
%           | ||x(t)||_1 - ||x(t-1)||_1 | / ||x(t)||_1 < rel_obj,
%       where x(t) is the estimate of the solution at iteration t.
%
%   - verbose: 0 no log, 1 a summary at convergence, 2 print main
%   steps (default: 1)
%
%   - weights: weights for a weighted L1-norm (default = 1)
%
%   - min_x/max_x: lower/upper bounds on x (default = {0, +Inf})
%
%   - mask_app: mask for the known fixed subpart of the image
%
%   - xA: fixed subpart of the image (default = 0)
%
%
%
% References:
% [1] M.J. Fadili and J-L. Starck, "Monotone operator splitting for
% optimization problems in sparse recovery" , IEEE ICIP, Cairo,
% Egypt, 2009.
% [2] Amir Beck and Marc Teboulle, "A Fast Iterative Shrinkage-Thresholding
% Algorithm for Linear Inverse Problems",  SIAM Journal on Imaging Sciences
% 2 (2009), no. 1, 183--202.
%
%

% Optional input arguments
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'Psit'), param.Psi = @(x) x; param.Psit = @(x) x; end
if ~isfield(param, 'tight'), param.tight = 1; end
if ~isfield(param, 'nu'), param.nu = 1; end
if ~isfield(param, 'rel_obj'), param.rel_obj = 1e-4; end
if ~isfield(param, 'max_iter'), param.max_iter = 200; end
if ~isfield(param, 'Psit'), param.Psit = @(x) x; end
if ~isfield(param, 'Psi'), param.Psi = @(x) x; end
if ~isfield(param, 'weights'), param.weights = 1; end
if ~isfield(param, 'pos'), param.pos = 0; end
if ~isfield(param, 'real'), param.real = 0; end
if ~isfield(param, 'min_x'), param.min_x = 0; end
if ~isfield(param, 'max_x'), param.max_x = +Inf; end
if ~isfield(param, 'mask_app'), param.mask_app = []; end
if ~isfield(param, 'xA'), param.xA = 0; end

psit_xA = param.Psit(param.xA) ;

% Useful functions
soft_ = @(z, T) sign(z).*max(abs(z)-T, 0);
soft = @(z,T) -psit_xA + soft_(z+psit_xA,T) ;

% Projection
if param.tight && ~param.pos && ~param.real % TIGHT FRAME CASE

temp = param.Psit(x);
sol = x + 1/param.nu * param.Psi(soft(temp, ...
                                      lambda*param.nu*param.weights)-temp);
crit_L1 = 'REL_OBJ'; iter_L1 = 1;
dummy = param.Psit(sol);
%coef=soft(temp,lambda*param.nu*param.weights);
norm_l1 = sum(param.weights(:).*abs(dummy(:)));

else % NON TIGHT FRAME CASE OR CONSTRAINT INVOLVED

% Initializations

sol = x;
% figure
% subplot 211
% imagesc(real(sol)), axis image, colorbar
% subplot 212
% imagesc(imag(sol)), axis image, colorbar
% pause
if param.pos || param.real
sol = real(sol);
end
% figure
% subplot 211
% imagesc(real(sol)), axis image, colorbar
% subplot 212
% imagesc(imag(sol)), axis image, colorbar
% pause
dummy = param.Psit(sol);
u_l1 = zeros(size(dummy));

prev_obj = 0; iter_L1 = 0;

% Soft-thresholding
% Init
if param.verbose > 1
fprintf('  Proximal L1 operator:\n');
end
while 1


% L1 norm of the estimate
norm_l1 = sum(param.weights(:).*abs(dummy(:)));

obj = .5*norm(x(:) - sol(:), 2)^2 + lambda * norm_l1;
rel_obj = abs(obj-prev_obj)/obj;

% Log
if param.verbose>1
fprintf('   Iter %i, prox_fval = %e, rel_fval = %e\n', ...
        iter_L1, obj, rel_obj);
end

% Stopping criterion
if (rel_obj < param.rel_obj)
crit_L1 = 'REL_OBJ'; break;
elseif iter_L1 >= param.max_iter
crit_L1 = 'MAX_IT_L1'; break;
end

% Soft-thresholding
res = u_l1*param.nu + dummy;
dummy = soft(res, lambda*param.nu*param.weights);
u_l1 = 1/param.nu * (res - dummy);
sol = x - param.Psi(u_l1);

if param.pos
sol = real(sol); 
sol(sol<param.min_x) = param.min_x(sol<param.min_x); 
sol_mask = sol(param.mask_app) ;
sol_mask(sol_mask>param.max_x) = param.max_x(sol_mask>param.max_x) ;
sol(param.mask_app) = sol_mask ;
end

if param.real
sol = real(sol);
end

% Update
prev_obj = obj;
iter_L1 = iter_L1 + 1;
dummy = param.Psit(sol);

end
end

% Log after the projection onto the L2-ball
if param.verbose >= 1
fprintf(['  prox_L1: prox_fval = %e,', ...
         ' %s, iter = %i\n'], norm_l1, crit_L1, iter_L1);
end

end
function [x, time_img] = imaging_fb_test(y, x, G0, sp_scale, Ax, N, K, nIter, param_algo)
%%

%% Image initialization from pre-calibrated DIEs

if ~isfield(param_algo, 'stop_crit'), param_algo.stop_crit = 1e-5; end
A = @(x) G0*so_fft2(x, K, sp_scale);
At = @(x) real(so_fft2_adj(G0'*x, N, K, sp_scale));
eta_l1 = param_algo.eta;
fprintf('================================================================\n');
fprintf('Imaging with true DDEs (eta = %5e)\n', eta_l1);
fprintf('================================================================\n');

% if strcmp(param_algo.opt_dict, 'Dirac')
%     prox_func =@(x) max(proxl1(x, eta_l1*(1.9/Ax) ), 0) ;
% else
%     param_l1.Psi = param_algo.Psi ;
%     param_l1.Psit = param_algo.Psit ;
%     param_l1.real = 1 ;
%     param_l1.pos = 1 ;
%     param_l1.verbose = 0 ;
%     param_l1.weights = 1 ;
%     prox_func =@(x) solver_prox_L1(x, eta_l1*(1.9/Ax), param_l1) ;
% end

param_l1.Psi = param_algo.Psi;
param_l1.Psit = param_algo.Psit;
param_l1.real = 1;
param_l1.pos = 1;
param_l1.verbose = 0;
param_l1.weights = 1;
prox_func = @(x) solver_prox_L1(x, eta_l1*(1.9/Ax), param_l1);

tic
for q = 1:nIter
    xold = x;
    grad = 2*At(A(x) - y);
    x = prox_func(x - (1.9/Ax)*grad);   
    crit_stop = norm(x - xold)/norm(xold);    
    if crit_stop < param_algo.stop_crit
        disp(['Stopping criterion reached: q = ', num2str(q)])
        break
    end
end
time_img = toc;

end

function [x, snr, crit] = imaging_fb(y, x_true, D1, D2, N, K, W, scale, u, v, S, J, nIter, param_algo)
%%

%% Image initialization from pre-calibrated DIEs

if ~isfield(param_algo, 'stop_crit'), param_algo.stop_crit = 1e-5; end
T = size(D1, 3);
G = cell(T, 1);
for t = 1:T
    G{t} = createGnufft_T2(D1(:, :, t), D2(:, :, t), [v(:, t), u(:, t)], K, S, J, W{t}); % see other option if needed...
end
G0 = cell2mat(G);
% st = nufft_init([v(:), u(:)], N, [J, J], K, N/2);
% G0 = st.p;
% scale = st.sn;

A = @(x) G0*so_fft2(x, K, scale);
At = @(x) real(so_fft2_adj(G0'*x, N, K, scale));
Ax = pow_method(A, At, size(x_true)); % operator norm, the factor 2 comes from the complex case
x = At(y)/Ax; % to be verified
eta_l1 = param_algo.eta;

fprintf('================================================================\n');
fprintf('Imaging with approximate DDEs (eta = %5e)\n', eta_l1);
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

param_l1.Psi = param_algo.Psi ;
param_l1.Psit = param_algo.Psit ;
param_l1.real = 1 ;
param_l1.pos = 1 ;
param_l1.verbose = 0 ;
param_l1.weights = 1 ;
prox_func = @(x) solver_prox_L1(x, eta_l1*(1.9/Ax), param_l1) ;

tic
for q = 1:nIter
    xold = x;
    grad = At(A(x) - y);
    x = prox_func(x - (1.9/Ax)*grad);   
    crit_stop = norm(x - xold)/norm(xold);
    if crit_stop < param_algo.stop_crit
        disp('Init. image from pre-calibrated DIEs')
        disp(['Stopping criterion reached: q = ', num2str(q)])
        disp(['norm(x-xold)/norm(x) = ', num2str(crit_stop)])
        break
    end
end
time_img = toc;

% Ouptut metrics
SNR = @(x, xtrue) 10*log10(sum(xtrue(:).^2)/sum((x(:) - xtrue(:)).^2));
crit = sum(abs(A(x) - y).^2) + eta_l1*sum(abs(param_algo.Psit(x)));
% error_dirty = sqrt(sum(sum((At(A(xinit) - y)).^2)));
snr = SNR(x, x_true);

% Display results
fprintf('%10s\t%15s\t%15s\t%15s\n', 'Iterations', 'Objective', 'SNR', 'Time')
fprintf('----------------------------------------------------------------\n');
fprintf('%10i\t%15e\t%15e\t%15e\n', q, crit, snr, time_img);
fprintf('================================================================\n');

end

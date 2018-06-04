%-------------------------------------------------------------------------%
% Main file to obtain the reference imaging results (with the true DDEs, 
% and normalized DIEs).
%-------------------------------------------------------------------------%
%% 
% [15/03/2018], P.-A. Thouvenin.
%-------------------------------------------------------------------------%

clc; clear all; close all;

addpath data
addpath utils
addpath lib
addpath lib/nufft
addpath src

im_choice = 'M31';
cov_type = 'vlaa';
test_number = '1';
test_aux = '';

rng(3);

% Data loading
name = [cov_type, '_', im_choice, '_', test_number];
name2 = [name, test_aux];
load(['synth_data_reg_eusipco_', name, '.mat'])

%% Imaging with the true DDEs

% -------------------------------------------------------------------------
% regularisation (wavelet dictionary)
% -------------------------------------------------------------------------
param_algo.nlevel = 4;
param_algo.opt_dict = 'SARA';
[param_algo.Psi, param_algo.Psit] = SARA_sparse_operator(randn(N), param_algo.nlevel,param_algo.opt_dict);
% -------------------------------------------------------------------------

% define algorithm parameters
x0 = zeros(size(x_th));
% % % param_algo.eta = [1e0, 5e-1, 1e-1, 5e-2, 1e-2, 5e-3, 1e-3, 5e-4, 1e-4, 5e-5, 1e-5, 5e-6, 1e-6]*2 ;
eta_range = [100, 50, 20, 10, 5, 1, 0.5];
nIter = 2e3; % 5e3

% Create G operator with the true DDEs
G = cell(T, 1);
parfor t = 1:T
    G{t} = createGnufft_T2(D_th(:, :, t), D_th(:, :, t), [v(:, t), u(:, t)], K, S, J, W{t}); % see other option if needed...
end
G0 = cell2mat(G);

A = @(x) G0*so_fft2(x, K, sp_scale);
At = @(x) real(so_fft2_adj(G0'*x, N, K, sp_scale));
Ax = 2*pow_method(A, At, size(x_th)); % operator norm, the factor 2 comes from the complex case
snr_ref = -Inf;
SNR = @(x, xtrue) 10*log10(sum(xtrue(:).^2)/sum((x(:) - xtrue(:)).^2));

% Run FB with mutliple values of eta, keep the version with the best SNR
for q = 1:length(eta_range)
    param_algo.eta = eta_range(q);
    [xq, snr_q, objective_q] = imaging_fb_th(y, x_th, G0, sp_scale, Ax, N, K, nIter, param_algo); % see time before convergence...
    x_true_dde = xq/max(xq(:));
    scale = sum(x_true_dde(:).*x_th(:))/sum(x_true_dde(:).^2);
    snr_true_dde = SNR(scale*x_true_dde, x_th);    
    if snr_true_dde >= snr_ref
        x = xq;
        snr_ref = snr_true_dde;
        snr = snr_q;
        objective = objective_q;
        eta = eta_range(q);
    end
end

x_true_dde = x/max(x(:));
scale = sum(x_true_dde(:).*x_th(:))/sum(x_true_dde(:).^2);
snr_true_dde = SNR(scale*x_true_dde, x_th);
disp(['SNR(x_true_dde) = ', num2str(snr_true_dde)]);

% results with approximate DDEs
nIter_FB = 2e3; % 1e4
D_approx = ones(1, na, T);
param_algo.eta = eta;
x_approx = imaging_fb(y, x_th, D_approx, D_approx, N, K, W, sp_scale, u, v, 1, J, nIter_FB, param_algo);

x_approx_die = x_approx/max(x_approx(:));
scale = sum(x_approx_die(:).*x_th(:))/sum(x_approx_die(:).^2);
snr_approx_die = SNR(scale*x_approx_die, x_th);
disp(['SNR(x_approx_die) = ', num2str(snr_approx_die)]);

mkdir('results')
save(['results/imaging_true_dde_', name2, '.mat'], 'x', 'snr', 'objective', 'eta', ...
      'x_true_dde', 'snr_true_dde', 'x_approx', 'x_approx_die', 'snr_approx_die', '-v7.3');

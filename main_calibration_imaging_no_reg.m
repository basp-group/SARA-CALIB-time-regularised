%-------------------------------------------------------------------------%
% Main file to perform joint calibration and imaging (results obtained from 
% main_imaging_true_dde.m are needed).
%-------------------------------------------------------------------------%
%% 
% [15/03/2018], P.-A. Thouvenin.
%-------------------------------------------------------------------------%

clc; clear all; close all;

addpath data
addpath utils;
addpath lib
addpath lib/nufft
addpath src
addpath results

im_choice = 'M31'; % M31, W28_256
cov_type = 'vlaa'; % vlaa, meerkat
test_number = '1';
test_aux = 'nu10';

rng(3);

% Load data and initialization variables
name = [cov_type, '_', im_choice, '_', test_number];
load(['synth_data_reg_eusipco_cell_', name, '.mat'])
load(['imaging_true_dde_', name, '.mat'], 'eta', 'x_approx')
name2 = [name, test_aux];
util_create_pool_bis(10, []); % 10
S2 = S*S;

% Regularisation parameters
param_algo.eta = eta;  % 1e-2 l1 image
param_algo.Jeps = 200; % 1000
param_algo.tol_x = 1e-5;
param_algo.tol_crit = 2e-2;

% Parameters for DDE estimation
param_dde.max_it = 50; % 200;
param_dde.JUtot = 5;   % 10;
param_dde.JU1 = 5;
param_dde.JU2 = 5;
param_dde.nu = 10;
mu = 0;
param_dde.tol_norm = 1e-5;

% Regularisation (wavelet dictionary)
param_algo.nlevel = 4;
param_algo.opt_dict = 'SARA';
[param_algo.Psi, param_algo.Psit] = SARA_sparse_operator(randn(N), param_algo.nlevel,param_algo.opt_dict);

%% Joint imaging and calibration (DDEs)

SNR = @(x, xtrue) 10*log10(sum(xtrue(:).^2)/sum((x(:) - xtrue(:)).^2));
x_approx_die = x_approx/max(x_approx(:));
scale = sum(x_approx_die(:).*x_th(:))/sum(x_approx_die(:).^2);
snr_approx_die = SNR(scale*x_approx_die, x_th);
disp(['SNR(x_approx_die) = ', num2str(snr_approx_die)]);

% DDE parameters
c = floor(S2/2) + 1;
param_dde.theta_maxR = 2*A*ones(S2,1); % 2*dde_var
param_dde.theta_minR = -param_dde.theta_maxR;
param_dde.theta_maxI = param_dde.theta_maxR;
param_dde.theta_minI = param_dde.theta_minR;
param_dde.theta_maxR(c) = (1 + 2*A_center); % 1 + 2*die_var
param_dde.theta_maxI(c) = 2*A_center;       % 2*die_var
param_dde.theta_minR(c) = (1 - 2*A_center); % 1 - 2*die_var
param_dde.theta_minI(c) = -2*A_center;

% Initialize parameters for image
x0 = x_approx;
x0(x0<0.05*max(x0(:))) = 0;
x0([1:5, end-5:end], :) = 0;
x0(:, [1:5, end-5:end]) = 0;
param_im.min_x = zeros(size(x0));
param_im.min_x(x0>0) = -5e-1*x0(x0>0);
param_im.max_x = 5e-1*x0(x0>0);
epsilon = zeros(size(x0));

% DDE initialisation: /!\ beware amplitude !
p = floor(P/2) + 1;
U = (A*(randn(S2, na, P) + 1i*randn(S2, na, P))/P)*sqrt(F);
U(:, :, p) = 0;
U(c, :, p) = sqrt(F);
D = computeD(U, F, [], [], T);
D1r = max(min(real(D), param_dde.theta_maxR), param_dde.theta_minR);
D1i = max(min(imag(D), param_dde.theta_maxI), param_dde.theta_minI);
D = D1r + 1i*D1i;

D1 = cell(T, 1);
for t = 1:T
    D1{t} = D(:, :, t);    
end
D2 = D1;

na_t = na*ones(T, 1);

% DDE estimation
[D1, D2, epsilon, objective, time_dde, error_dirty, snr_dirty, snr_x] = joint_calibration_imaging(y, Y, x0, ...
    epsilon, x_th, u, v, Omega, na_t, T, J, K, N, S, D1, D2, sp_scale, ...
    V, param_im, param_dde, param_algo, W, name2, mu);

x_no_reg = x0 + epsilon;
x_no_reg = x_no_reg/max(x_no_reg(:));
scale = sum(x_no_reg(:).*x_th(:))/sum(x_no_reg(:).^2);
snr_no_reg = SNR(scale*x_no_reg, x_th);
disp(['SNR(x_no_reg) = ', num2str(snr_no_reg)]);
% figure; plot(objective(1:find(objective, 1, 'last')));

%% Save the results

mkdir('results')
save(['results/joint_imaging_dde_no_reg_', name2, '.mat'], 'D1', 'D2', 'epsilon', ...
      'x0', 'objective', 'param_dde', 'time_dde', 'error_dirty', 'snr_dirty', 'snr_x', ...
      'snr_no_reg', 'x_no_reg', 'snr_approx_die', 'x_approx_die', '-v7.3');

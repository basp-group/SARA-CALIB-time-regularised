%-------------------------------------------------------------------------%
% Main file to perform joint calibration and imaging with temporal 
% regularization (results obtained from main_imaging_true_dde.m are needed).
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
addpath results

im_choice = 'M31'; % W28_256 M31
cov_type = 'vlaa'; % vlaa, meerkat
test_number = '1';
test_aux = 'nu10';

rng(3);

% Loading data and initialization variables
name = [cov_type, '_', im_choice, '_', test_number];
load(['synth_data_reg_eusipco_', name,'.mat'])
load(['imaging_true_dde_', name, '.mat'], 'eta', 'x_approx')
name2 = [name, test_aux];
util_create_pool_bis(10, []);
S2 = S*S;

% Regularisation parameters
param_algo.eta = eta;
param_algo.Jeps = 200; % 1000
param_algo.tol_x = 1e-5;
param_algo.tol_crit = 2e-2;

% Parameters for DDE estimation
param_dde.max_it = 50; % 200;
param_dde.JUtot = 5; % 10;
param_dde.JU1 = 5;   % 5
param_dde.JU2 = 5;   % 5
param_dde.nu = 10; % 1000;
mu = 0;
param_dde.tol_norm = 1e-5;

% regularisation (wavelet dictionary)
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
p = floor(P/2) + 1;
theta_maxR = 2*A*ones(S2,1);
theta_minR = -theta_maxR;
theta_maxI = theta_maxR;
theta_minI = theta_minR;
theta_maxR(c) = (1 + 2*A_center);  
theta_maxI(c) = 2*A_center;         
theta_minR(c) = (1 - 2*A_center);   
theta_minI(c) = -2*A_center;

param_dde.theta_maxR = 2*A*ones(S2, P);   
param_dde.theta_minR = -2*A*ones(S2, P);
param_dde.theta_maxR(c, p) = 1 + 2*A;
param_dde.theta_minR(c, p) = 1 - 2*A;
param_dde.theta_maxI = 2*A*ones(S2, P);
param_dde.theta_minI = -param_dde.theta_maxI;

% Initialize image parameters
x0 = x_approx;
x0(x0<0.05*max(x0(:))) = 0 ;
x0([1:5, end-5:end], :) = 0;
x0(:, [1:5, end-5:end]) = 0; 
param_im.min_x = zeros(size(x0));
param_im.min_x(x0>0) = -5e-1*x0(x0>0);
param_im.max_x = 5e-1*x0(x0>0);
epsilon = zeros(size(x0));

% DDE initialisation: /!\ beware amplitude
p = floor(P/2) + 1;
U1 = (A*(randn(S2, na, P) + 1i*randn(S2, na, P))/P);
U1(:, :, p)=  0;
U1(c, :, p) = 1;
U1r = max(min(real(U1), reshape(param_dde.theta_maxR, [S2, 1, P])), reshape(param_dde.theta_minR, [S2, 1, P]));
U1i = max(min(imag(U1), reshape(param_dde.theta_maxI, [S2, 1, P])), reshape(param_dde.theta_minI, [S2, 1, P]));
U1 = U1r + 1i*U1i;
U2 = U1;

% DDE estimation (set parameters & estimate DDE + image)
[U1, U2, epsilon, objective, time_dde, error_dirty, snr_dirty, snr_x] = joint_calibration_imaging_reg(y, Y, x0, ...
    epsilon, x_th, u, v, Omega, na, T, J, P, F, K, N, S, U1, ...
    U2, sp_scale, V, param_im, param_dde, param_algo, W, name2, mu);

x_reg = x0 + epsilon;
x_reg = x_reg/max(x_reg(:));
scale = sum(x_reg(:).*x_th(:))/sum(x_reg(:).^2);
snr_reg = SNR(scale*x_reg, x_th);
disp(['SNR(x_reg) = ', num2str(snr_reg)]);
% figure; plot(objective(1:find(objective, 1, 'last')));

%% Save the results

mkdir('results')
save(['results/joint_imaging_dde_reg_', name2,'.mat'], 'U1', 'U2', 'epsilon', ...
     'x0', 'objective', 'param_dde', 'time_dde', 'error_dirty', 'snr_dirty', 'snr_x', ...
     'snr_reg', 'x_reg', 'snr_approx_die', 'x_approx_die', '-v7.3');

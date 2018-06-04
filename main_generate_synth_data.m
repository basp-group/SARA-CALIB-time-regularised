%-------------------------------------------------------------------------%
% Main file to generate the synthetic datasets (different formats for the 
% input data, depending on the presence or absence of temporal 
% in the solver).
%-------------------------------------------------------------------------%
%% 
% [15/03/2018], P.-A. Thouvenin.
%-------------------------------------------------------------------------%

clc; clear all; close all;
format compact;

addpath utils 
addpath lib
addpath lib/nufft 
addpath src
addpath synth

rng(3)
util_create_pool_bis(4, []);

% Image parameters / name
im_choice = 'M31'; % possible choices: M31, W28, W28_256, W28_512
cov_type = 'vlaa'; % possible choices: 'vlaa', 'meerkat'
test_number = '1';
N = 256*[1,1];     % image sizes
K = 2*N;           % size of the spatial Fourier domain
name = [cov_type, '_', im_choice, '_', test_number];
dataName_cell = ['data/synth_data_reg_eusipco_cell_', name]; % name for the cell format dataset
dataName = ['data/synth_data_reg_eusipco_', name];

% Coverage parameters
T = 200;
na = 27;
M = T*na*(na-1)/2;
n_pairs = na*(na-1)/2;
dl = 1.2; % pixel size (typically 1.5)
input_snr = 150;  % noise level 

% DDE parameters
S = 5; % size of the DDE spatial support (S^2 elements in total)
J = 5; % size of the gridding kernels (J^2 elements in total)
P = 3; % size of the DDE temporal support (3, ..., 11)
F = T; % size of the temporal Fourier space
A = 4e-2; % off-center DDE amplitude
A_center = 5e-2; % DDE amplitude (center elements in the temporal and spatial Fourier spaces)

% Data generation
[Y_, V_, W, Omega_, u_ab, v_ab, y, x_th, sp_scale, U_, D_, N, K, na] = generate_synth_data(S, J, N, K, P, T, F, dl, input_snr, cov_type, A, im_choice);
mkdir('data')

%% Save data for the standard calibration algorithm (modify data format)
u = cell(T, 1);
v = cell(T, 1);
V = cell(T, 1);
D_th = cell(T, 1);
Omega = cell(T, 1);
y = y(:);
Y = cell(T, 1);
for t = 1:T
    Y{t} = Y_(:, :, t);
    Omega{t} = Omega_(:, :, :, t);
    u{t} = u_ab(:, t);
    v{t} = v_ab(:, t);
    V{t} = V_(:, :, :, t);
    D_th{t} = D_(:, :, t);
    U_th = U_;
end
save(dataName_cell, 'x_th', 'y', 'Y', 'Omega', 'V', 'sp_scale', 'D_th', 'U_th', 'u', 'v', 'T', 'na', 'N', 'K', 'S', 'J', 'P', 'F', 'W', 'A', 'A_center');

%% Save data for the time-regularized calibration algorithm
u = u_ab;
v = v_ab;
V = V_;
D_th = D_;
U_th = U_;
Omega = Omega_;
Y = Y_;
y = y(:);
save(dataName, 'x_th', 'y', 'Y', 'Omega', 'V', 'sp_scale', 'D_th', 'U_th', 'u', 'v', 'T', 'na', 'N', 'K', 'S', 'J', 'P', 'F', 'W', 'A', 'A_center');

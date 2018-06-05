function [U1, U2, epsilon, objective, time_tot, error_dirty, snr_dirty, snr_x] = joint_calibration_imaging_reg(y, Y, x0, epsilon, x_th, u, v, Omega, na, T, J, P, F, K, N, S, U1, U2, scale, V, param_im, param_dde, param_algo, W, name, mu)
% Joint imaging and calibration procedure in absence of temporal
% regularization.
%-------------------------------------------------------------------------%
%%
% Input:
% > y       : data vector [M*T, 1], M = na*(na-1)/2, na: number of
%             antennas, T: number of snapshots 
% > Y       : data with redundancy [na, na, T]
%             (Nan on the diag., y_ab = conj(y_ab) for a > b)
% > x0      : reference image
% > epsilon : initial perturbation (generally zeros(size(x0)))
% > x_th    : ground truth image
% > u       : u component [M, T]
% > v       : v component [M, T]
% > Omega   : u-v components with redundancy (similar to Y) [na, na, 2, T],
%             Omega(:, :, 1, t) -> v(t)
%             Omega(:, :, 2, t) -> u(t)
% > na      : number of antennas
% > T       : number of snapshots
% > J       : size of the gridding kernels (square kernels [J, J])
% > P       : size of the DDE support in the temporal Fourier space
% > F       : size of the temporal Fourier space (F = T)
% > K       : size of the spatial Fourier space [2, 1]
% > N       : image dimension [2, 1]
% > S       : spatial dimension of the DDE kernels (square kernels [S, S])
% > U1/U2   : initial DDE kernels (in the temporal Fourier domain) [S2, na, P]
% > scale   : scaling coefficients from the NUFFT (size(scale) = N)
% > V       : spatial NUFFT gridding kernels with redundancy (similar to Y) [J2, na, na, T]
% > param_im  : image-related parameters ()
% > param_dde : DDE-related parameters (contraints definition)
%               - .theta_maxR (.theta_maxI): maximum value for the real
%                                             (imag.) part
%               - .theta_minR (.theta_minI): minimum value for the real
%                                             (imag.) part
%               - .max_it :maximum number of iteration of the malor loop
%               - .JUtot : number of iteration for the DDEs
%               - .JU1 : number of iteration to estimate U1
%               - .JU2 : number of iteration to estimate U2
%               - .nu : hyperparameter controlling the Eucl. distance
%                        between U1 and U2
%               - .tol_norm : stopping criterion of the DDE major loop
% > param_algo : algorithm parameters
%                - .min_x / .max_x : bounds on the perturbation espilon
% > W          : spatial NUFFT gridding kernels [J2, T]
% > name       : name of the temporary backup file
% > mu         : additional hyperparameter (see eq (13) in Audrey's SPIE paper)
%
% Output:
% < U1/U2       : estimated DDEs [S2, na, P]   
% < epsilon     : estimated pertubration
% < objective   : objective function
% < time_tot    : computation time for each iteration of the outermost loop
% < error_dirty : SNR of the dirty image
% < snr_dirty   : SNR of the dirty image
% < snr_x       : image reconstruction quality (in dB)
%
%-------------------------------------------------------------------------%
%% 
% [15/03/2018], P.-A. Thouvenin.
%-------------------------------------------------------------------------%
%% Remarks:
% - add code to implement temporal NUFFT (if applicable)

%% Initialization
time_tot = zeros(param_dde.max_it+1, 1);
objective = zeros(param_dde.max_it+1, 1);
error_dirty = zeros(param_dde.max_it+1, 1);
snr_dirty = zeros(param_dde.max_it+1, 1);
snr_x = zeros(param_dde.max_it+1, 1);

param_l1.min_x = param_im.min_x;
param_l1.max_x = param_im.max_x;
param_l1.mask_app = (x0 > 0);
param_l1.Psi = param_algo.Psi;
param_l1.Psit = param_algo.Psit;
param_l1.real = 1;
param_l1.pos = 1;
param_l1.verbose = 0;
param_l1.weights = 1;
param_l1.xA = x0;

lambda_scale = 1/prod(K);
param_algo.eta = param_algo.eta*lambda_scale;
proxEpsilon = @(x,Ax) solver_prox_L1_full_image(x, param_algo.eta*1.9/Ax, param_l1);
regul_x = @(x) param_algo.eta*sum( abs(param_l1.Psit(x0+x)) ) ;
SNR = @(x, xtrue) 10*log10(sum(xtrue(:).^2)/sum((x(:) - xtrue(:)).^2));

theta_maxoR = param_dde.theta_maxR;
theta_minoR = param_dde.theta_minR;
theta_maxoI = param_dde.theta_maxI;
theta_minoI = param_dde.theta_minI;
JU2o = param_dde.JU2;
JU1o = param_dde.JU1;
nuo = param_dde.nu;

% keyboard
    
% Create Gt operator
Gt = []; % unsued parameter (left for legacy purpose, to be cleansed later...)
scale_t = [];
S2 = S^2;

% Initial value of the objective function
G = cell(T, 1);
D1 = computeD(U1, T, Gt, scale_t, T);
D2 = computeD(U2, T, Gt, scale_t, T);

parfor t = 1:T
    G{t} = createGnufft_T2(D1(:, :, t), D2(:, :, t), [v(:, t), u(:, t)], K, S, J, W{t}); % incorporate misisng antenna positions (evolves w.r.t. t)
end
G0 = cell2mat(G);
clear G

B_tmp = @(x) G0*so_fft2(x, K, scale);
Bt_tmp = @(x) real(so_fft2_adj(G0'*x, N, K, scale));
Phi = zeros(size(U1));
Phi(floor(S2/2)+1, :, floor(P/2)+1) = 1;

data_fid = lambda_scale*sum(abs(B_tmp(x0 + epsilon) - y).^2);
err1 = nuo*sum(abs(U1(:) - U2(:)).^2)/2;
err2 = mu*(sum(abs(U1(:) - Phi(:)).^2) + sum(abs(U2(:) - Phi(:)).^2))/2; 
objective(1) = data_fid(1) + err1 + err2 + regul_x(epsilon);
error_dirty(1) = sqrt(sum(sum((Bt_tmp(B_tmp(x0 + epsilon) - y)).^2)));
snr_dirty(1) = SNR(Bt_tmp(B_tmp(x0 + epsilon)), Bt_tmp(y));
snr_x(1) = SNR(x0 + epsilon, x_th);

%% Global algorithm
fprintf('===============================================================================================\n');
fprintf('Joint calibration and imaging (temporal regularization) [%s]\n', name);
fprintf('===============================================================================================\n');
tic
for nIter = 1:param_dde.max_it
    %%% Update DDEs
    x_hat = reshape(fftshift(reshape(so_fft2(x0 + epsilon, K, scale), K)), [prod(K), 1]);
    
    % Computation of the portions of x_hat to be transmitted to the workers
    id = true(na, 1);
    X1 = cell(na, 1);
    X2 = cell(na, 1);
    for a = 1:na % relatively time consuming -> apparently no viable alternative for the later data distribution, no parfor allowed here (memory considerations)
        id(a) = false;
        idx1 = [];
        idx2 = [];
        for t = 1:T
            % for U1
            id_a = find(~isnan(squeeze(Y(a, id, t))));
            om1 = squeeze(Omega(a, id, 1, t)).'; % /!\ squeeze does not remove the first singleton dimension in this case -> need to transpose)
            om1 = om1(id_a);
            om2 = squeeze(Omega(a, id, 2, t)).';
            om2 = om2(id_a);
            idt = unique(indices4Xhat(S, J, K, [om1, om2])); % see if I can directly generate the indices of interest
            idx1 = unique([idx1; idt]); % [S2^2*J2*(na-1), na]
            % for U2
            id_a = find(~isnan(squeeze(Y(id, a, t))));
            om1 = squeeze(Omega(id, a, 1, t)); % check dimensions
            om1 = om1(id_a);
            om2 = squeeze(Omega(id, a, 2, t));
            om2 = om2(id_a);
            idt = unique(indices4Xhat(S, J, K, [om1, om2])); 
            idx2 = unique([idx2; idt]);
        end
        id(a) = true;
        X1{a} = sparse(idx1, ones(size(idx1)), x_hat(idx1), prod(K), 1); % fine, no specfic overhead when reading entries from a sparse matrix
        X2{a} = sparse(idx2, ones(size(idx2)), x_hat(idx2), prod(K), 1);
    end
    clear om1 om2 x_hat
    
    for nIter_dde = 1:param_dde.JUtot
        fprintf('--DDE iterations: %4i / %4i \n', nIter_dde, param_dde.JUtot)
        
        % Update U1 & U2
        U_old = U1;
        parfor a = 1:na % see parfeval for the update of U1 and U2: compute indices, then send appropriate portion of x_hat
            % update U1a
            id = true(na, 1);
            id(a) = false;  
            Ya = squeeze(Y(a, :, :));
            Va = squeeze(V(:, a, :, :));
            Omega_a = squeeze(Omega(a,:,:,:)); % [na, 1, 2, T] om1 = squeeze(Omega(a,id,1,t)).';
            % build H1_a [na-1, S2, T], n = na-1
            Ha = zeros(na-1, S2, T);
            for t = 1:T % most time consuming operation
                id_a = find(~isnan(Ya(id, t))); % find Y(a, id, t)
                %idx = indices4xhat_2(ll0, kxJ, kyJ, K); % generate indices from the bloc start/end positions previously obtained 
                % [S2^2*J2*(na-1), na] % /!\ nonzeros systematically returns a column vector
                om1 = squeeze(Omega_a(id,1,t));
                om1 = om1(id_a);
                om2 = squeeze(Omega_a(id,2,t));
                om2 = om2(id_a);
                idt = indices4Xhat(S, J, K, [om1, om2]); % [S2^2*J2*(na-1), 1]
                if ~isempty(id_a)
                    D_at = D2(:, id, t);       % [S2, na-1, 1]
                    Ha(id_a, :, t) = createH1a(full(X1{a}(idt)), Va(:, id, t), D_at, J, S, id_a).'; % [n, S2] 
                end
            end
            U1(:, a, :) = updateU1a(Ya(id, :), reshape(U1(:, a, :), [S2, P]), ...
                          reshape(U2(:, a, :), [S2, P]), Ha, Gt, scale_t, F, T, JU1o, nuo, mu, ...
                          lambda_scale, theta_minoR, theta_maxoR, theta_minoI, theta_maxoI);
        end

        % update D1
        D1 = computeD(U1, T, Gt, scale_t, T);
        stop_dde = norm(U1(:) - U_old(:))/norm(U_old(:));
        
        % update U2a
        U_old = U2;
        parfor a = 1:na   
            id = true(na, 1);
            id(a) = false; 
            Ya = squeeze(Y(:, a, :));
            Va = squeeze(V(:, :, a, :));
            Omega_a = squeeze(Omega(:,a,:,:)); % [na, 1, 2, T]
            % build H2_a, [na-1, S2, T], n = na-1
            Ha = zeros(na-1, S2, T);
            for t = 1:T
                id_a = find(~isnan(Ya(id, t)));
                om1 = squeeze(Omega_a(id,1,t));
                om1 = om1(id_a);
                om2 = squeeze(Omega_a(id,2,t));
                om2 = om2(id_a);
                idt = indices4Xhat(S, J, K, [om1, om2]); % [S2^2*J2*(na-1), 1] % /!\ nonzeros systematically returns a column vector
                if ~isempty(id_a)
                    D_at = D1(:, id, t); % [S2, na-1, 1]
                    Ha(id_a, :, t) = createH2a(full(X2{a}(idt)), squeeze(Va(:, id, t)), D_at, J, S, id_a); % -> revoir format d'entr�e + Va
                end
            end
            U2(:, a, :) = updateU2a(Ya(id, :), reshape(U1(:, a, :), [S2, P]), ...
                          reshape(U2(:, a, :), [S2, P]), Ha, Gt, scale_t, F, T, JU2o, nuo, mu, ...
                          lambda_scale, theta_minoR, theta_maxoR, theta_minoI, theta_maxoI);
        end
        % update D2
        D2 = computeD(U2, T, Gt, scale_t, T);
       
        % DDE stopping criterion
        stop_dde = max(stop_dde, norm(U2(:) - U_old(:))/norm(U_old(:)));
        stop_crit = max(stop_dde) ;
        if stop_crit(end) < param_dde.tol_norm
            fprintf('--DDE iterations: stopping criterion satisfied');
            break
        end        
    end 
    fprintf('-----------------------------------------------------------------------------------------------\n');
    clear X1 X2
    
    %%% Update Image
    
    % Create convolution matrix G
    G = cell(T, 1);
    parfor t = 1:T % group several time instants if needed...
        G{t} = createGnufft_T2(D1(:, :, t), D2(:, :, t), [v(:, t), u(:, t)], K, S, J, W{t});
    end
    G0 = cell2mat(G);
    
    % Update epsilon
    B_tmp = @(x) G0*so_fft2(x, K, scale);
    Bt_tmp = @(x) real(so_fft2_adj(G0'*x, N, K, scale));
    epsilon = updateEpsilon(y, x0, epsilon, param_algo.Jeps, param_algo.tol_x, nIter, proxEpsilon, B_tmp, Bt_tmp, lambda_scale); % [28/02/2018] OK, there was a problem in the solver I had (wrong versions on the WS, who modified this?!)
    
    % Display results
    time_tot(nIter+1) = toc;
    err1 = nuo*sum(abs(U1(:) - U2(:)).^2)/2;
    err2 = mu*(sum(abs(U1(:) - Phi(:)).^2) + sum(abs(U2(:) - Phi(:)).^2))/2; 
    data_fid = lambda_scale*sum(abs(B_tmp(x0 + epsilon) - y).^2);
    objective(nIter+1) = data_fid + err1 + err2 + regul_x(epsilon);
    error_dirty(nIter+1) = sqrt(sum(sum((Bt_tmp(B_tmp(x0 + epsilon) - y)).^2)));
    snr_dirty(nIter+1) = SNR(Bt_tmp(B_tmp(x0 + epsilon)), Bt_tmp(y));
    snr_x(nIter+1) = SNR(x0 + epsilon, x_th);
    
    x_reg = x0 + epsilon;
    x_reg = x_reg/max(x_reg(:));
    scale_x_reg = sum(x_reg(:).*x_th(:))/sum(x_reg(:).^2);
    scaled_snr = SNR(scale_x_reg*x_reg, x_th);
    
    fprintf('%10s\t%15s\t%15s\t%15s\t%15s\t%15s\n', 'Iteration', 'Objective', 'SNR', 'scaled SNR', 'SNR (dirty)', 'Time')
    fprintf('-----------------------------------------------------------------------------------------------\n');
    fprintf('%10i\t%15e\t%15e\t%15e\t%15e\t%15e\n', nIter, objective(nIter+1), snr_x(nIter+1), ...
             scaled_snr, snr_dirty(nIter+1), time_tot(nIter+1));
    fprintf('===============================================================================================\n');
    
    if mod(nIter, 5)==0
        save(['results/img_dde_reg', name], 'U1', 'U2', 'epsilon')
    end
    
    % Global stopping criterion
    if abs((objective(nIter+1) - objective(nIter))/objective(nIter)) < param_algo.tol_crit
        fprintf('End of algorithm: global stopping criterion satisfied.\n')
        break
    end
end

end

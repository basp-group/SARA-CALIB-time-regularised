function [U1, U2, epsilon, objective, time_tot, error_dirty, snr_dirty, snr_x] = joint_calibration_imaging(y, Y, x0, epsilon, x_th, u, v, Omega, na, T, J, K, N, S, U1, U2, scale, V, param_im, param_dde, param_algo, W, name, mu)
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
% > x       : ground truth image
% > u       : u component [M, T]
% > v       : v component [M, T]
% > Omega   : u-v components with redundancy (similar to Y) [na, na, 2, T],
%             Omega(:, :, 1, t) -> v(t)
%             Omega(:, :, 2, t) -> u(t)
% > na      : number of antennas
% > T       : number of snapshots
% > J       : size of the gridding kernels (square kernels [J, J])
% > K       : size of the spatial Fourier space [2, 1]
% > N       : image dimension [2, 1]
% > S       : spatial dimension of the DDE kernels (square kernels [S, S])
% > U1/U2   : initial DDE kernels {T, 1}, [S2, na] for each cell
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

% < U1/U2       : estimated DDEs {T, 1}, [S2, na] for each cell    
% < epsilon     : estimated pertubration
% < objective   : objective function
% < time_tot    : computation time for each iteration of the outermost loop
% < error_dirty : SNR of the dirty image
% < snr_dirty   : SNR of the dirty image
% < snr_x       : image reconstruction quality (in dB)
%-------------------------------------------------------------------------%
%% 
% [15/03/2018], P.-A. Thouvenin.
%-------------------------------------------------------------------------%

%% Notes:
% J: size of the NFFT support
% S: size of the DDE support
% N: image dimension
% K: Fourier space dimension (spatial)
% Omega: {T,1}: cell containing the u, v coordinates with redundancy (as for Y matrices, Omega{t}: [na(t), na(t), 2], first slice corresponding to v)
% V: cell containing the gridding coefficients for NFFT, presents redundancy as for Y
% scale: scaling factors associated to the NFFT
% U1 / U2: cell containing the nonzero coefficients of the DDE

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
proxEpsilon = @(x,Ax) solver_prox_L1_full_image(x, param_algo.eta*1.9/Ax, param_l1) ;
regul_x = @(x) param_algo.eta * sum( abs(param_l1.Psit(x0+x)) ) ;
SNR = @(x, xtrue) 10*log10(sum(xtrue(:).^2)/sum((x(:) - xtrue(:)).^2));

theta_maxoR = param_dde.theta_maxR;
theta_minoR = param_dde.theta_minR;
theta_maxoI = param_dde.theta_maxI;
theta_minoI = param_dde.theta_minI;
JU2o = param_dde.JU2;
JU1o = param_dde.JU1;
nuo = param_dde.nu;
stop_crit_temp = zeros(T,1);
err_temp = zeros(T,1);
S2 = S^2;
Phi = sparse((floor(S2/2)+1)*ones(na(1), 1), (1:na(1)).', ones(na(1), 1), S2, na(1));

% Initial value of the objective function
G = cell(T, 1);
parfor t = 1:T % parfor
%     id_upper = getIndicesUp(na(t)); % indices corresponding to the upper half of the square matrix V{t}
    G{t} = createGnufft_T2(U1{t}, U2{t}, [v{t}, u{t}], K, S, J, W{t}); % V{t}(:, id_upper) see st.uu..., see format
    err_temp(t) = nuo*sum(abs(U1{t}(:) - U2{t}(:)).^2) + mu*(sum(abs(U1{t}(:) - Phi(:)).^2) + sum(abs(U2{t}(:) - Phi(:)).^2));
end
G0 = cell2mat(G);
B_tmp = @(x) G0*so_fft2(x, K, scale);
Bt_tmp = @(x) real(so_fft2_adj(G0'*x, N, K, scale));

data_fid = lambda_scale*sum(abs(B_tmp(x0 + epsilon) - y).^2);
objective(1) = data_fid + sum(err_temp)/2 + regul_x(epsilon);
error_dirty(1) = sqrt(sum(sum((Bt_tmp(B_tmp(x0 + epsilon) - y)).^2)));
snr_dirty(1) = SNR(Bt_tmp(B_tmp(x0 + epsilon)), Bt_tmp(y));
snr_x(1) = SNR(x0 + epsilon, x_th);

fprintf('===============================================================================================\n');
fprintf('Joint calibration and imaging (no temporal regularization) [%s]\n', name);
fprintf('===============================================================================================\n');
%% Global algorithm
tic
for nIter = 1:param_dde.max_it
    
    %%% Update DDEs
    x_hat = reshape(fftshift(reshape(so_fft2(x0 + epsilon, K, scale), K)), [prod(K), 1]);
    
    % Selection of the portions of x_hat to be transmitted to the workers
    X1 = cell(T, 1);
    X2 = cell(T, 1);
    for t = 1:T % relatively time consuming -> apparently no viable alternative for the later data distribution, no parfor allowed here (memory considerations)
        id = true(na(t),1);
        idx1 = [];
        idx2 = [];
        for a = 1:na(t)
            id(a) = false;
            % for U1
            id_a = find(~isnan(squeeze(Y{t}(a, id))));
            om1 = squeeze(Omega{t}(a, id, 1)).'; % /!\ squeeze does not remove the first singleton dimension in this case -> need to transpose)
            om1 = om1(id_a);
            om2 = squeeze(Omega{t}(a, id, 2)).';
            om2 = om2(id_a);
            idt = unique(indices4Xhat(S, J, K, [om1, om2])); % see if I can directly generate the indices of interest
            idx1 = unique([idx1; idt]); % [S2^2*J2*(na-1), na]
            % for U2
            id_a = find(~isnan(squeeze(Y{t}(id, a))));
            om1 = squeeze(Omega{t}(id, a, 1)); % check dimensions
            om1 = om1(id_a);
            om2 = squeeze(Omega{t}(id, a, 2));
            om2 = om2(id_a);
            idt = unique(indices4Xhat(S, J, K, [om1, om2])); 
            idx2 = unique([idx2; idt]);
            id(a) = true;
        end
        X1{t} = sparse(idx1, ones(size(idx1)), x_hat(idx1), prod(K), 1); % fine, no specfic overhead when reading entries from a sparse matrix
        X2{t} = sparse(idx2, ones(size(idx2)), x_hat(idx2), prod(K), 1);
    end
    clear om1 om2 x_hat
    
    for nIter_dde = 1:param_dde.JUtot
        fprintf('--DDE iterations: %4i / %4i \n', nIter_dde, param_dde.JUtot);
        % Update U1 & U2      
        parfor t = 1:T % see parfeval or parfor for the update of U1 and U2: compute indices, then send appropriate portion of x_hat
            U_old = U1{t};          
            id = true(na(t),1);
            Ht = zeros(na(t)-1, S2, na(t));
            for a = 1:na(t)
                id(a) = false;
                id_nnz = find(~isnan(Y{t}(a,id)));
                om1 = squeeze(Omega{t}(a, id, 1)).'; % /!\ squeeze does not remove the first singleton dimension in this case -> need to transpose)
                om1 = om1(id_a);
                om2 = squeeze(Omega{t}(a, id, 2)).';
                om2 = om2(id_a);
                idt = indices4Xhat(S, J, K, [om1, om2]); % see if I can directly generate the indices of interest
                Ht(:, :, a) = createH1a(full(X1{t}(idt)), squeeze(V{t}(:, a, id)), U2{t}(:,id), J, S, id_nnz).';         
                id(a) = true;
            end 
            U1{t} = updateU1t(Y{t}, U1{t}, U2{t}, Ht, na(t), JU1o, nuo, mu, lambda_scale, theta_minoR, theta_maxoR, theta_minoI, theta_maxoI);
            stop_crit_temp(t) = norm(U1{t} - U_old, 'fro')/norm(U_old, 'fro');
            
            % update U2{t}
            U_old = U2{t};
            id = true(na(t),1);
            Ht = zeros(na(t)-1, S2, na(t));
            for a = 1:na(t)
                id(a) = false;
                id_nnz = find(~isnan(Y{t}(id,a)));
                om1 = squeeze(Omega{t}(id, a, 1)); % check dimensions
                om1 = om1(id_a);
                om2 = squeeze(Omega{t}(id, a, 2));
                om2 = om2(id_a);
                idt = indices4Xhat(S, J, K, [om1, om2]); 
                Ht(:, :, a) = createH2a(full(X2{t}(idt)), squeeze(V{t}(:, id, a)), U1{t}(:,id), J, S, id_nnz);  % [S2^2*J2*(na-1), na]
                id(a) = true;
            end 
            U2{t} = updateU2t(Y{t}, U1{t}, U2{t}, Ht, na(t), JU2o, nuo, mu, lambda_scale, theta_minoR, theta_maxoR, theta_minoI, theta_maxoI);
            stop_crit_temp(t) = max(stop_crit_temp(t), norm(U2{t} - U_old, 'fro')/norm(U_old, 'fro'));
            
            % regularization term
            err_temp(t) = nuo*sum(abs(U1{t}(:) - U2{t}(:)).^2) + mu*(sum(abs(U1{t}(:) - Phi(:)).^2) + sum(abs(U2{t}(:) - Phi(:)).^2));
        end

        % DDE stopping criterion
        stop_crit = max(stop_crit_temp) ;
        if stop_crit(end) < param_dde.tol_norm
            fprintf('--DDE iterations: stopping criterion satisfied');
            break
        end
    end 
    fprintf('-----------------------------------------------------------------------------------------------\n');
    
    %% Update Image
    
    % Create convolution matrix G
    G = cell(T,1);
    parfor t = 1:T
        % id_upper = getIndicesUp(na(t)); % indices corresponding to the upper half of the square matrix V{t}
        G{t} = createGnufft_T2(U1{t}, U2{t}, [v{t}, u{t}], K, S, J, W{t}); % V{t}(:,id_upper)
    end
    G0 = sparse(cell2mat(G));

    % Update epsilon
    B_tmp = @(x) G0*so_fft2(x, K, scale);
    Bt_tmp = @(x) real(so_fft2_adj(G0'*x, N, K, scale));
    epsilon = updateEpsilon(y, x0, epsilon, param_algo.Jeps, param_algo.tol_x, nIter, proxEpsilon, B_tmp, Bt_tmp, lambda_scale);
    
    % Display monitoring results
    time_tot(nIter+1) = toc;
    data_fid = lambda_scale*sum(abs(B_tmp(x0 + epsilon) - y).^2);
    objective(nIter+1) = data_fid + sum(err_temp)/2 + regul_x(epsilon);
    error_dirty(nIter+1) = sqrt(sum(sum((Bt_tmp(B_tmp(x0 + epsilon) - y)).^2)));
    snr_dirty(nIter+1) = SNR(Bt_tmp(B_tmp(x0 + epsilon)), Bt_tmp(y));
    snr_x(nIter+1) = SNR(x0 + epsilon, x_th);
    
    x_reg = x0 + epsilon;
    x_reg = x_reg/max(x_reg(:));
    scale_x = sum(x_reg(:).*x_th(:))/sum(x_reg(:).^2);
    scaled_snr = SNR(scale_x*x_reg, x_th);
    
    fprintf('%10s\t%15s\t%15s\t%15s\t%15s\t%15s\n', 'Iteration', 'Objective', 'SNR', 'scaled SNR', 'SNR (dirty)', 'Time')
    fprintf('-----------------------------------------------------------------------------------------------\n');
    fprintf('%10i\t%15e\t%15e\t%15e\t%15e\t%15e\n', nIter, objective(nIter+1), snr_x(nIter+1), ...
             scaled_snr, snr_dirty(nIter+1), time_tot(nIter+1));
    fprintf('===============================================================================================\n');
   
    if mod(nIter, 5)==0 % to be desactivated if needed
        save(['results/img_dde_no_reg_', name], 'U1', 'U2', 'epsilon')
    end
    
    % Global stopping criterion
    if abs((objective(nIter+1) - objective(nIter))/objective(nIter) ) < param_algo.tol_crit
        fprintf('End of algorithm: global stopping criterion satisfied.\n')
        break
    end
end

end

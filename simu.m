%%%%%%%%%% Simulation %%%%%%%%%%%%
% Model: yi(s) = g(xi.T*beta(s)) + zi.T*alpha(s) + eta_i(s) + epsilon_i(s)

clear all;
addpath('./utilities/');
addpath('/SIVC_2016_code/')

%% Set Parameters
nv = 200;
n = 100;
nsimu = 100;
% p = 3;
% G = [1,1/2,1; 0,1,4];
G0 = [1, 1/2, -1];
ul = 3;
index_g = 1; % index_g = 2;
seed = 2021;

sigma_xi = [0.5 0.5]/5;
sigma_epsilon = 0.05;

% generate data 
[beta0, alpha0, alpha_star, S, simu_x, simu_y, simu_eta, ...
    simu_eta_star, simu_g, simu_z] = simu_data_gen(n, nv, nismu, ... 
        G0, ul, index_g, sigma_xi, sigma_epsilon, seed);


%% Estimation errors for each simulation dataset
% PLSIM
all_betaest = zeros(p,nv,nsimu);
all_alphaest = zeros(p,nv,nsimu);
all_lambdaest = zeros(nv,nsimu);
all_gest = zeros(n,nv,nsimu);
all_dgest = zeros(n,nv,nsimu);
all_etaest = zeros(n,nv,nsimu);
all_betaini =zeros(p,nv,nsimu);
% competitors: SIVC & MVCM
all_betaest_SIVC = zeros(p,nv,nsimu);
all_betaest_MVCM = zeros(p,nv,nsimu);
all_gest_SIVC = zeros(n,nv,nsimu);
all_dgest_SIVC = zeros(n,nv,nsimu);

% settings in estimation
smooth_coeff = 1;
verbose = 0; 

tic;
for nn = 1:nsimu
% can use parallel computing to speed it up 
% for nn = 1
    rng(nn);
    x = simu_x(:,:,nn);
    ally = simu_y(:,:,nn);

    disp('Estimate coefficients and index function');
    [betaest, alphaest, lambdaest, g_est, dg_est, etaest, beta_ini] = EST(x, ally,...
        S, smooth_coeff, verbose, nn);
    toc;
    all_betaest(:,:,nn) = betaest;
    all_alphaest(:,:,nn) = alphaest;
    all_lambdaest(:,nn) = lambdaest;
    all_gest(:,:,nn) = g_est;
    all_dgest(:,:,nn) = dg_est;
    all_etaest(:,:,nn) = etaest;
    all_betaini(:,:,nn) = beta_ini;

    disp('SIVC estimators');
    for s=1:nv
        all_betaest_SIVC(:, s, nn) = getBeta(x, ally, nv, S, s, hx, hg, h, beta_ini(:, s));
    end
    allxb = x*all_betaest_SIVC(:, :, nn); %X.T * beta, (n,m)
    for i=1:n
        for s=1:nv
            xb0 = x(i, :)*all_betaest_SIVC(:, s, nn); 
            [all_gest_SIVC(i, s, nn), all_dgest_SIVC(i,s,nn)] = locallinear1(1, xb0, hg, allxb(:), ally(:));
        end
    end
    toc;    

    disp('MVCM estimators');
    for s=1:nv
        options = optimoptions('fminunc', 'Algorithm','Quasi-Newton', 'Display', 'off');
        B0 = zeros(p, 2);
        B0(:, 1) = all_betaest_SIVC(:, s, nn); B0(:, 2) = 0;
        Best = fminunc(@(B) wls(x, ally, S', B, h2, s), B0, options);
        all_betaest_MVCM(:, s, nn) = Best(:, 1);
    end
    toc;

end 


%% Mean Estimation Errors 
% (i) PLSIM
% beta(s) 
tmp = abs(all_betaest - repmat(beta0, 1, 1, nsimu)); tmp = mean(tmp, 3);
MIAE_beta = mean(tmp, 2) % mean integrated absolute error for beta
tmp = (all_betaest - repmat(beta0, 1, 1, nsimu)).^2; tmp = mean(tmp, 3); %tmp = tmp.^0.5;
MISE_beta = mean(tmp, 2) % mean integrated squared error for beta
% alpha(s)
tmp = abs(all_alphaest - repmat(alpha_star, 1, 1, nsimu)); tmp = mean(tmp, 3);
MIAE_alpha = mean(tmp, 2) % mean integrated absolute error for beta
tmp = (all_alphaest - repmat(alpha_star, 1, 1, nsimu)).^2; tmp = mean(tmp, 3); %tmp = tmp.^0.5;
MISE_alpha = mean(tmp, 2) % mean integrated squared error for beta
% g
tmp = abs(all_gest - simu_g); tmp = mean(tmp, 2);  tmp = mean(tmp, 3);
MIAE_g = mean(tmp, 1) % mean integrated absolute error for g
tmp = (all_gest - simu_g).^2; tmp = mean(tmp, 2);  tmp = mean(tmp, 3);%tmp = tmp.^0.5;
MISE_g = mean(tmp, 1) % mean integrated squared error for g

% (ii) SIVC 
% beta(s)
tmp = abs(all_betaest_SIVC - repmat(beta0, 1, 1, nsimu)); tmp = mean(tmp, 3);
MIAE_beta_SIVC = mean(tmp, 2) % mean integrated absolute error for beta
tmp = (all_betaest_SIVC - repmat(beta0, 1, 1, nsimu)).^2; tmp = mean(tmp, 3); %tmp = tmp.^0.5;
MISE_beta_SIVC = mean(tmp, 2)
% g
tmp = abs(all_gest_SIVC - simu_g); tmp = mean(tmp, 2);  tmp = mean(tmp, 3);
MIAE_g_SIVC = mean(tmp, 1) % mean integrated absolute error for g
tmp = (all_gest_SIVC - simu_g).^2; tmp = mean(tmp, 2);  tmp = mean(tmp, 3);%tmp = tmp.^0.5;
MISE_g_SIVC = mean(tmp, 1)

% (iii) MVCM
tmp = abs(all_betaest_MVCM - repmat(beta0, 1, 1, nsimu)); tmp = mean(tmp, 3);
MIAE_beta_MVCM = mean(tmp, 2) % mean integrated absolute error for beta
tmp = (all_betaest_MVCM - repmat(beta0, 1, 1, nsimu)).^2; tmp = mean(tmp, 3); %tmp = tmp.^0.5;
MISE_beta_MVCM = mean(tmp, 2)


%% save the estimators & errors 
estname = sprintf('est_result_n%d_s%d_nsimu%d_g%d_ul%d.mat', n,nv,nsimu,index_g,ul)
save(estname, 'all_betaest', 'all_alphaest', 'all_lambdaest', 'all_gest', ...
    'all_dgest', 'all_etaest', 'all_betaini', 'all_betaest_SIVC', 'all_gest_SIVC',...
    'all_dgest_SIVC', 'all_betaest_MVCM', 'smooth_coeff', 'MIAE_beta', ...
    'MISE_beta', 'MIAE_alpha', 'MISE_alpha', 'MIAE_g', 'MISE_g', 'MIAE_beta_SIVC', ...
    'MISE_beta_SIVC', 'MIAE_g_SIVC', 'MISE_g_SIVC', 'MIAE_beta_MVCM', 'MISE_beta_MVCM');


%% Prediction Errors 
disp('prediction Errors');
pps = [0.3, 0.5, 0.7]; % the proportion for dividing the data
all_SE = zeros(3,nsimu, length(pps));
all_AE = zeros(3,nsimu, length(pps));
tic;
for idx = 1:length(pps)
    pp = pps(idx);
    fprintf('proportion=%.1f\n', pp)
    for nn = 1:nsimu
        nn 
        x = simu_x(:,:,nn);
        ally = simu_y(:,:,nn);
        
        [SE, AE] = prediction(x, ally, S, smooth_coeff, verbose, pp, nn);
        all_SE(:,nn,idx) = SE;
        all_AE(:,nn,idx) = AE;
        toc;
    end 
end 
mean_SE = squeeze(mean(all_SE,2))
mean_AE = squeeze(mean(all_AE,2))
save(estname, 'all_SE', 'all_AE', 'mean_SE', 'mean_AE', 'pps', '-append');



%% Identifying the space of hidden factors Z
% (i) residual matrix & its projection
all_R = simu_y - all_gest;
all_U = zeros(n,n,nsimu);
for nn = 1:nsimu
    R = all_R(:,:,nn);
    [U, ~, ~] = svd(R);
    all_U(:,:,nn) = U;

    x = simu_x(:,:,nn);
    Projx = x * (x'*x)^(-1) * x'; % (n,n)
    PR = Projx * R; % (n,nv)
    all_Projx(:,:,nn) = Projx;
    all_PR(:,:,nn) = PR;
end
resname = sprintf('residuals_n%d_s%d_nsimu%d_g%d_ul%d.mat', n,nv,nsimu,index_g,ul)
save(resname, 'all_R', 'all_U', 'all_Projx', 'all_PR', 'q');


% (ii) factor analysis to obtain q_est 
% use R package for factor analysis to determine which method to use
% then obtian q_est by the chosen method
q_est = 1;


% (iii) correlation between U_{1:q} & Z
pval_Z = zeros(q_est,q_est,nsimu);
corr_UZ = zeros(q_est,q_est,nsimu);    
for nn = 1:nsimu
    U = all_U(:,:,nn);
    [corr_UZ(:,:,nn), pval_Z(:,:,nn)] = corr(U(:,1:q_est), simu_z(:,:,nn));
end 
corr_UZ = squeeze(corr_UZ);
pval_Z = squeeze(pval_Z);
mean(corr_UZ)
std(corr_UZ)
save(resname, 'q_est'，'corr_UZ', 'pval_Z', '-append');

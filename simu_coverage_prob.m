%%%%%%%%%% Simulation: Coverage probs for SCB of beta(s) & g(x*beta(s)) %%%%%%%%%%%%
% Model: yi(s) = g(xi.T*beta(s)) + zi.T*alpha(s) + eta_i(s) + epsilon_i(s)
clear all;
addpath('./utilities/');

%% Set Parameters
nv = 200;
n = 100;
nsimu = 100;
G0 = [1, 1/2, -1];
ul = 3;
index_g = 1; % index_g = 2;
seed = 2021;
sigma_xi = [0.5 0.5]/5;
sigma_epsilon = 0.05;

SCB_alpha = [0.05, 0.01]; % significant levels

%% generate data & load estimators 
[beta0, ~, ~, S, simu_x, simu_y, ~, ~, simu_g, ~] = simu_data_gen(n, ...
	nv, nismu, G0, ul, index_g, sigma_xi, sigma_epsilon, seed);

estname = sprintf('est_result_n%d_s%d_nsimu%d_g%d_ul%d.mat', n,nv,nsimu,index_g,ul)
load(estname, 'all_betaest', 'all_alphaest', 'all_lambdaest', 'all_gest',...
	'all_etaest', 'smooth_coeff')


%% Bootstrap for SCB
disp('Obtain simultaneous confidence bands');
verbose = 0; 
B = 500;

all_Cb = zeros(p, nsimu, length(SCB_alpha));
all_Cg = zeros(nsimu, length(SCB_alpha));
tic;
for nn = 1:nsimu
    nn 
    x = simu_x(:,:,nn);
    ally = simu_y(:,:,nn);
    alphaest = all_alphaest(:,:,nn);
    betaest = all_betaest(:,:,nn);
    lambdaest = all_lambdaest(:,nn);
    g_est = all_gest(:,:,nn);
    etaest = all_etaest(:,:,nn);

    [Cb, Cg, ~, ~] = SCB(x, ally, S, alphaest, betaest, lambdaest, g_est, ...
        etaest, smooth_coeff, verbose, SCB_alpha, B, nn);
    all_Cb(:,nn,:) = Cb;
    all_Cg(nn,:) = Cg;
    toc;
end 


%% Summarize the bounds for each replication
disp('Obtain the lower bound & upper bound')
all_betaL = zeros(p,nv,2,nsimu);
all_betaU = zeros(p,nv,2,nsimu);
all_sort_g0 = zeros(n*nv, nsimu);
all_sort_gest = zeros(n*nv, nsimu);
all_sort_gL = zeros(n*nv, 2, nsimu);
all_sort_gU = zeros(n*nv, 2, nsimu);
for nn = 1:nsimu
    % SCB for beta(s)
    all_betaL(:,:,:,nn) = repmat(all_betaest(:,:,nn), 1,1,2) - repmat(all_Cb(:,nn,:), 1,nv,1); % (p,nv,2)
    all_betaU(:,:,:,nn) = repmat(all_betaest(:,:,nn), 1,1,2) + repmat(all_Cb(:,nn,:), 1,nv,1); % (p,nv,2)
    
    % SCB for g
    x = simu_x(:,:,nn); 
    betaest = all_betaest(:,:,nn);
    xb = x * betaest; 
    [sort_xb, index] = sort(xb(:));
    sort_g0 = simu_g(:,:,nn); 
    sort_g0 = sort_g0(index); % (n*nv, 1)
    all_sort_g0(:,nn) = sort_g0; 
    % all_sort_g0(:,:,nn) = g(sort_xb); % (n*nv, 1)
    sort_gest = all_gest(:,:,nn);
    sort_gest = sort_gest(index); % (n*nv, 1)
    all_sort_gest(:,nn) = sort_gest; 
    all_sort_gL(:,:,nn) = repmat(sort_gest, 1,2) - all_Cg(nn,:); % (n*nv, 2)
    all_sort_gU(:,:,nn) = repmat(sort_gest, 1,2) + all_Cg(nn,:); % (n*nv, 2)
end 


%% calculate the coverage probs 
disp('Calculate the coverage prob')
idx = 1:nsimu;
idx0 = idx(all_Cg(:,1)~=0);

diff_bL = repmat(beta0, 1,1,2,nsimu) - all_betaL; % (p,nv,2,nsimu)
diff_bU = all_betaU - repmat(beta0, 1,1,2,nsimu); % (p,nv,2,nsimu)
% cover_b = mean(diff_bL >= 0, 2) .* mean(diff_bU >= 0, 2); % (p,1,2,nsimu)
cover_b = mean(diff_bL(:,:,:,idx0) >= 0, 2) .* mean(diff_bU(:,:,:,idx0) >= 0, 2); % (p,1,2,nsimu)
cover_b = mean(squeeze(cover_b)==1, 3) % (p, 2)

diff_gL = repmat(all_sort_g0, 1,2,1) - all_sort_gL; % (n*nv, 2, nsimu)
diff_gU = all_sort_gU - repmat(all_sort_g0, 1,2,1); % (n*nv, 2, nsimu)
% cover_g = mean(diff_gL >= 0, 1) .* mean(diff_gU >= 0, 1); % (1, 2, nsimu)
mid_idx = (0.02 * n * nv + 1) : ((1-0.02) * n * nv);
cover_g = mean(diff_gL(mid_idx,:,idx0) >= 0, 1) .* mean(diff_gU(mid_idx,:,idx0) >= 0, 1); % (1, 2, nsimu)
cover_g = mean(squeeze(cover_g)==1,2) % (2,1)


%% save the results
resname = sprintf('cover_prob_n%d_s%d_nsimu%d_g%d_ul%d.mat', n,nv,nsimu,index_g,ul)
save(resname, 'all_Cg', 'all_Cb', 'B', 'SCB_alpha', 'cover_g', 'cover_b');

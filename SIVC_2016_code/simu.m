clear all;

%% Path
% pwd = 'E:\SDR_test_rev03092010';
% pwd = '/Users/xinchaoluo/Academic/Softwares/Matlab/Functions';
% addpath(genpath(pwd));

%% Set parameters
%rng(19880531)

n = 200;
n = 40;
% p = 4;
p = 2;
m = 50;
rho = 0.6;
lambda = [1 0.5];
% sigma = 0.3;
nsimu = 200; % number of simulation
G = 500;

muX = zeros(1, p);
SigmaX = zeros(p, p);
for i=1:p
    for j=1:p
        SigmaX(i, j) = rho^(abs(i-j));
    end
end

beta = zeros(p, m);
t = linspace(0, 1, m);
for s=1:m
    beta(1, s) = 1 + t(s)^2;
    beta(2, s) = (1 - t(s))^2;
    % beta(3, s) = 4*t(s)*(1 - t(s));
    % beta(4, s) = -1 + 4*(t(s)-0.5)^2;
    beta(:, s) = beta(:, s)/sqrt(sum(beta(:, s).^2));
end

hx = n^(-1/3)*rho^0.5;   %cn^(-1/3)
hy = n^(-1/5)*rho^0.5;   %cn^(-1/5)
h = m^(-1/5)*0.3;
% h1 = 0.3;
% h2 = 0.15;
%% Generate data
simu_ally = zeros(n, m, nsimu);
simu_x = zeros(n, p, nsimu);
simu_eta = zeros(n, m, nsimu);
simu_g = zeros(n, m, nsimu);

% g & eta functions
eta1 = @(s) sqrt(2)*sin(2*pi.*s); %psi1(s)
eta2 = @(s) sqrt(2)*cos(2*pi.*s); %psi2(s)
g = @(xb) sin(2*xb) + 2*cos(2 + xb); %g1(xb)

for nn=1:nsimu % for nn_th simulation
    x =  mvnrnd(muX, SigmaX, n);
    simu_x(:, :, nn) = x;
    % eta(s) = sum_{k=1}^{2}[xi_jk*psi_k(s)]
    simu_eta(:, :, nn) = [sqrt(lambda(1))*randn(n,1),sqrt(lambda(2))*randn(n,1)]*[eta1(t') eta2(t')]';
    for s=1:m  
        xb = x*beta(:, s); % X.T * beta(s), (n,1)
        simu_g(:, s, nn) = g(xb); % g1(xb), (n,1)
        simu_ally(:, s, nn) = simu_g(:, s, nn) + simu_eta(:, s, nn); % + sigma*randn(n, 1); %y(s) = g(xb)+eta(s)
    end
end

%% Estimation

all_betaest = zeros(p, m, nsimu);
all_gest = zeros(n, m, nsimu);
all_Sigma_etaest = zeros(m, m, nsimu); %Simga_eta, (m,m)
% all_sigma2est = zeros(1, nsimu);
tic
for nn=1:nsimu
% for nn=1:2
    nn
    x = simu_x(:, :, nn);
    ally = simu_ally(:, :, nn);  
    
% inital beta    
%     beta0 = zeros(p, m);  
%     for s=1:m
%         y=ally(:, s);
%         [WX, W, f, d] = ldr(y, x, 'LAD', 'cont', 1, 'nslices', 5);
%         W = W*sign(W(1));
%         beta0(:, s) = W;
%     end 
%  % smooth
%     tmp = zeros(p, m);
%     for j=1:p
%         for s=1:m
%             tmp(j, s) = locallinear0(1, t(s), 0.15, t', beta0(j, :)');
%         end        
%     end
%     
%     beta0 = tmp;

    beta0 = beta; % starting point
    
    
% Estimate beta
    for s=1:m
           s
        % beta(sm) at grid sm, (p,1)
        all_betaest(:, s, nn) = getBeta(x, ally, m, t, s, hx, hy, h, beta0(:, s));
%         all_betaest(:, s, nn) = all_betaest(:, s, nn)*sign(all_betaest(1, s, nn));
    end
%         all_betaest(:, :, nn) = SIVC_sif(all_betaest(:, :, nn), t);
%     for s=1:m
%         all_betaest(:, s, nn) = all_betaest(:, s, nn)/sqrt(sum(all_betaest(:, s, nn).^2));
%     end
    
% Estimate g
    allxb = x*all_betaest(:, :, nn); %X.T * beta, (n,m)
    h1 = cvh1(x, all_betaest(:, :, nn), allxb, ally); % compute h1 using CV
    for i=1:n
        for s=1:m
            xb0 = x(i, :)*all_betaest(:, s, nn); %xi.T * beta(sm), scalar
            % solve eq(13), view Zim=1 to estimate g_hat(X.T*beta(s))
            all_gest(i, s, nn) = locallinear0(1, xb0, h1, allxb(:), ally(:));
            % allxb(:),(n*m,1); ally(:),(n*m,1);
        end
    end
    
% Estimate eta & sigma
    ystar = ally - all_gest(:, :, nn); % epsilon_star = y-g(X.t*beta) = eta+error, (n,m)
    etaest = zeros(n, m);
    h2 = cvh2(t, ystar); % compute h2 using CV
    for i=1:n
        for s=1:m
            % solve eq(14), view Wms=1 to estimate eta_i(s)
            etaest(i, s) = locallinear0(1, t(s), h2, t', ystar(i, :)');        
        end
    end
    all_Sigma_etaest(:, :, nn) = etaest'*etaest/(n-p); % eq(16), compute covariance matrix
    % res = ally - all_gest(:, :, nn) - etaest; % residual, epsilon, (n,m)
    % all_sigma2est(nn) = sum(sum(res.^2))/n/m; % sigma^2 of epsilon
   
    toc;
end
toc

%% Inference
n_used = 2; % the n_th simulation
gest = all_gest(:, :, n_used);
betaest = all_betaest(:, :, n_used);
% compute the scalar C_beta, C_g for simultaneous confidence bands
% ystar = simu_ally(:,:,n_used) - gest;
% [varCb, varCg] = SCB(ystar, gest, x, m, t, G, betaest, hx, hy, h, h1); 

%% Cover Probobility
[biasbeta, biasg] = allBias(x, ally, ystar, betaest, m, t, hx, hy, h, h1); %compute the bias
% [cp1, cp2] = copr(beta, betaest, biasbeta, simu_g(:, :, n_used), gest, biasg, varCb, varCg);
%% Analysis
betaest = mean(all_betaest, 3);
Sigma_etaest = mean(all_Sigma_etaest, 3);
position = [0.05 0.7 0.17 0.2; 0.3 0.7 0.17 0.2; 0.05 0.35 0.17 0.2; 0.3 0.35 0.17 0.2];

% beta
figure;
for kk = 1:p
    subplot(2,p/2,kk);
    subplot('Position',position(kk, :));
    % use eq(17) to construct simultaneous confidence bands for beta
    % betaU = betaest(kk, :) - biasbeta(kk, :) + varCb(kk); % upper bound for beta
    % betaL = betaest(kk, :) - biasbeta(kk, :) - varCb(kk); % lower bound for beta
    plot(t, beta(kk,:)', 'k-'); hold on; % ground truth
    plot(t, betaest(kk,:)', 'k:', 'lineWidth', 1.4); % mean of #nsimu beta_est
%     plot(t, betaU', 'k-.', 'lineWidth', 1.4);
%     plot(t, betaL', 'k-.', 'lineWidth', 1.4);
    %axis([0 1 -1 2]);
    title(strcat('\beta_',num2str(kk)));
    xlabel('s');
end


% g
xb_simu = zeros(n, m, nsimu);
for nn=1:nsimu
    for s=1:m
        xb_simu(:, s, nn) = simu_x(:, :, nn)*all_betaest(:, s, nn);
    end
end
sort_xb = xb_simu(:); 
sort_g = g(sort_xb); 
[sort_xb, index] = sort(sort_xb);
sort_g = sort_g(index);
sort_gest = all_gest(:);
sort_gest = sort_gest(index);
figure;
plot(sort_xb, sort_g, 'k-'); hold on;
plot(sort_xb, sort_gest, 'k:', 'lineWidth', 1.4);

xb_simu = zeros(n, m);
for s=1:m
    xb_simu(:, s) = simu_x(:, :, n_used)*all_betaest(:, s, n_used); 
    % ground truth of X.T * beta(s), (n,1)
end
sort_xb = xb_simu(:); % vectorize X.T * beta, (n*m,1)
sort_g = g(sort_xb); % vectorized ground truth g(X.T*beta), (n*m,1)
[sort_xb, index] = sort(sort_xb);
sort_g = sort_g(index);
sort_gest = all_gest(:, :, n_used) - biasg; %(n,m)
sort_gest = sort_gest(:); % vectorize, (n*m,1)
sort_gest = sort_gest(index);

% figure;
subplot('Position',[0.55 0.35 0.4 0.55]);
% use eq(18) to construct simultaneous confidence bands for g
% gU = sort_gest + varCg; % upper bound for g
% gL = sort_gest - varCg; % lower bound for g
plot(sort_xb, sort_g, 'k-'); hold on; % ground truth
plot(sort_xb, sort_gest, 'k:', 'lineWidth', 1.4); % estimate for n_uesd-th simulation
% plot(sort_xb, gU, 'k-.', 'lineWidth', 1.4); 
% plot(sort_xb, gL, 'k-.', 'lineWidth', 1.4);
axis([-4 3 -4 5])
title('g');
xlabel(strcat('X^T\beta'));
   

% Sigma_eta
Sigma_eta = 2*lambda(1)*sin(2*pi*t)'*sin(2*pi*t) + 2*lambda(2)*cos(2*pi*t)'*cos(2*pi*t);
% ground truth of Sigma_eta using eq(16), (m,m)
rl = min(Sigma_eta(:)); ru = max(Sigma_eta(:));

figure;
subplot(2,2,1);
imagesc(Sigma_eta, [rl, ru]); colorbar;title('true \Sigma_\eta(s,t)')
subplot(2,2,2);
imagesc(Sigma_etaest, [rl, ru]); colorbar; title('estimated \Sigma_\eta(s,t)')


% lambda&phi
[phi, lambda1] = eig(Sigma_eta); % eigen decomposition of ground truth Sigma_eta
lambda1 = diag(lambda1);
[lambda1, index] = sort(lambda1, 'descend');
phi = phi(:, index); phi(:, 2) = -phi(:, 2);
phi(:,1) = -phi(:,1);

all_phiest = zeros(m, m, nsimu);
all_lambdaest = zeros(2, nsimu);
tmp1 = eta1(t(10))./phi(10, 1);
tmp2 = eta2(t(10))./phi(10, 2);
for nn=1:nsimu
    [phi2, lambda2] = eig(all_Sigma_etaest(:, :, nn));
    lambda2 = diag(lambda2);
    [lambda2, index] = sort(lambda2, 'descend');
    all_phiest(:, :, nn) = phi2(:, index); 
    all_phiest(:, 1, nn) = all_phiest(:, 1, nn)*sign(all_phiest(10, 1, nn)); 
    all_lambdaest(1, nn) = lambda2(1)/tmp1^2;
    all_phiest(:, 2, nn) = all_phiest(:, 2, nn)*sign(all_phiest(1, 2, nn)); 
    all_lambdaest(2, nn) = lambda2(2)/tmp2^2;
end
phi2 = mean(all_phiest, 3);

% plot psi_1(s)
subplot(2,2,3);
plot(phi(:, 1), 'k-'); hold on; 
plot(phi2(:, 1), 'k:', 'lineWidth', 1.2); 
title('\psi_1(s)');
% plot psi_2(s)
subplot(2,2,4);
plot(phi(:, 2), 'k-'); hold on; 
plot(phi2(:, 2), 'k:', 'lineWidth', 1.2); 
title('\psi_2(s)');


% ME
tmp = abs(all_betaest - repmat(beta, 1, 1, nsimu)); tmp = mean(tmp, 3);
MIAE_beta = mean(tmp, 2) % mean integrated absolute error for beta

tmp = (all_betaest - repmat(beta, 1, 1, nsimu)).^2; tmp = mean(tmp, 3); %tmp = tmp.^0.5;
MISE_beta = mean(tmp, 2) % mean integrated squared error for beta

tmp = abs(all_gest - simu_g); tmp = mean(tmp, 2);  tmp = mean(tmp, 3);
MIAE_g = mean(tmp, 1) % mean integrated absolute error for g

tmp = (all_gest - simu_g).^2; tmp = mean(tmp, 2);  tmp = mean(tmp, 3);%tmp = tmp.^0.5;
MISE_g = mean(tmp, 1) % mean integrated squared error for g

tmp = abs(all_phiest - repmat(phi, 1, 1, nsimu)); tmp = mean(tmp, 3);
MIAE_phi = mean(tmp, 1); MIAE_phi = MIAE_phi(1:2)

tmp = (all_phiest - repmat(phi, 1, 1, nsimu)).^2; tmp = mean(tmp, 3); %tmp = tmp.^0.5;
MISE_phi = mean(tmp, 1); MISE_phi = MISE_phi(1:2)

tmp = abs(all_lambdaest - repmat(lambda, nsimu, 1)');
MAE_lambda = mean(tmp, 2) % mean absolute error for lambda

tmp = (all_lambdaest - repmat(lambda, nsimu, 1)').^2;
RMSE_lambda = mean(tmp, 2).^0.5 % root mean square error for lambda

% tmp = abs(all_sigma2est - repmat(sigma^2, 1, nsimu));
% MAE_sigma2 = mean(tmp, 2)
% 
% tmp = (all_sigma2est - repmat(sigma^2, 1, nsimu)).^2;
% RMSE_sigma2 = mean(tmp, 2).^0.5

% save n200m50true10.50.3(19880531h0.3e-1).mat
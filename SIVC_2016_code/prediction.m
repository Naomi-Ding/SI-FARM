clear all;

%% Path
pwd = 'D:\Dropbox\SDR_test_rev03092010';
% pwd = '/Users/xinchaoluo/Academic/Softwares/Matlab/Functions';
% pwd = '/home/xluo/SVIM/SDR_test_rev03092010';
addpath(genpath(pwd));

%% Set parameters
rng(19880531)

n = 200;
p = 4;
m = 50;
rho = 0.6;
pp = 0.7;
lambda = [1 0.5];
sigma = 0.3;
nsimu = 200;


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
    beta(3, s) = 4*t(s)*(1 - t(s));
    beta(4, s) = -1 + 4*(t(s)-0.5)^2;
    beta(:, s) = beta(:, s)/sqrt(sum(beta(:, s).^2));
end

hx = n^(-1/3)*rho^0.5;
hy = n^(-1/5)*rho^0.5;
h = m^(-1/5)*0.3;
h1 = 0.3;
h2 = 0.15;

%% Generate data
simu_ally = zeros(n, m, nsimu);
simu_x = zeros(n, p, nsimu);
simu_eta = zeros(n, m, nsimu);
simu_g = zeros(n, m, nsimu);

% g & eta functions
eta1 = @(s) sqrt(2)*sin(2*pi.*s);
eta2 = @(s) sqrt(2)*cos(2*pi.*s);
g = @(xb) sin(2*xb) + 2*cos(2 + xb);

for nn=1:nsimu
    x =  mvnrnd(muX, SigmaX, n);
    simu_x(:, :, nn) = x;
    simu_eta(:, :, nn) = [sqrt(lambda(1))*randn(n,1),sqrt(lambda(2))*randn(n,1)]*[eta1(t') eta2(t')]';
    for s=1:m  
        xb = x*beta(:, s);
        simu_g(:, s, nn) = g(xb);
        simu_ally(:, s, nn) = simu_g(:, s, nn) + simu_eta(:, s, nn) + sigma*randn(n, 1);
    end
end

%% Estimation
all_betaest = zeros(p, m, nsimu);
all_gest = zeros(n, m, nsimu);
all_Sigma_etaest = zeros(m, m, nsimu);
all_sigma2est = zeros(1, nsimu);

SE = zeros(2, nsimu);
AE = zeros(2, nsimu);

tic
for nn = 1:nsimu
    nn
    x = simu_x(:, :, nn);
    ally = simu_ally(:, :, nn);  
    
    n_tra = floor(n*pp);
    n_tes = n - n_tra;
    ind = 1:n;
    ind_tra = randsample(n, n_tra);
    ind_tes = ind;
    ind_tes(ind_tra) = [];

    x_tra = x(ind_tra, :); ally_tra = ally(ind_tra, :);
    x_tes = x(ind_tes, :); ally_tes = ally(ind_tes, :);
    
    
%% Estimation

    
% inital beta    
    beta0 = zeros(p, m);  
    for s=1:m
        y=ally(:, s);
        [WX, W, f, d] = ldr(y, x, 'LAD', 'cont', 1, 'nslices', 5);
        W = W*sign(W(1));
        beta0(:, s) = W;
    end 
 % smooth
    tmp = zeros(p, m);
    for j=1:p
        for s=1:m
            tmp(j, s) = locallinear0(1, t(s), h2, t', beta0(j, :)');
        end        
    end
    
    beta0 = tmp;
%     beta0 = beta;
    
    %SVIM
    betaest_simu = zeros(p, m);
    for s=1:m
        betaest_simu(:, s) = getBeta(x_tra, ally_tra, m, t, s, hx, hy, h, beta0(:, s));
    end

    ally_tes_est = zeros(n_tes, m);
    allxb = x_tes*betaest_simu;
    for i=1:n_tes
        for s=1:m
            xb0 = x_tes(i, :)*betaest_simu(:, s);
            ally_tes_est(i, s) = locallinear0(1, xb0, h2, allxb(:), ally_tes(:));        
        end
    end
    SE(1, nn) = (sum(sum((ally_tes_est - ally_tes).^2))/n_tes/m)^0.5;
    AE(1, nn) = (sum(sum(abs(ally_tes_est - ally_tes)))/n_tes/m);
    
    %MVCM
    ally_tes_est2 = zeros(n_tes, m);
    betaest_simu2 = zeros(p, m);
%     betaest_tmp = zeros(p, m);
%     for s=1:m
%         betaest_tmp(:, s) = (x_tra'*x_tra)\(x_tra'*ally_tra(:, s));
%     end
%     for j=1:p
%         for s=1:m
%             delta = t - t(s);
%             kernel = kh(delta, h2);
%             betaest_simu2(j, s) = betaest_tmp(j, :)*kernel/sum(kernel);           
%         end
%     end

    for s=1:m
%         options = optimoptions('fsolve', 'Algorithm','levenberg-marquardt', 'Display', 'off');
%         options = optimoptions('fsolve', 'Display', 'off');
        options = optimoptions('fminunc', 'Algorithm','Quasi-Newton', 'Display', 'off');
        B0 = zeros(p, 2);
        B0(:, 1) = betaest_simu(:, s); B0(:, 2) = 0;
%         Best = fsolve(@(B) wls(x_tra, ally_tra, t, B, h2, s), B0, options);
        Best = fminunc(@(B) wls(x_tra, ally_tra, t, B, h2, s), B0, options);
        betaest_simu2(:, s) = Best(:, 1);
    end
    
    for s=1:m
        ally_tes_est2(:, s) = x_tes* betaest_simu2(:, s);
    end
    SE(2, nn) = (sum(sum((ally_tes_est2 - ally_tes).^2))/n_tes/m)^0.5;
    AE(2, nn) = (sum(sum(abs(ally_tes_est2 - ally_tes)))/n_tes/m);
end
toc

%% Analysis
mean(AE, 2)
std(AE, 0, 2)
mean(SE, 2)
std(SE, 0, 2)

% save simucompare0.7dr.mat

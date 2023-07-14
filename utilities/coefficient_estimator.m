function [betaest, alphaest, lambdaest, beta_ini] = coefficient_estimator(x, y, S, hx, hg, hs, lambda_ini, smooth_coeff, verbose)

%% Input: 
    % S: (nv, 1) vector, the vector of grids for the measurement of response
    % lambda_ini: (1,nv), random vecotr from U(0,1)
    % smooth_coeff: boolean, whether smooth the estimators (alpha(s), beta(s))
    % verbose: whether displaying the convergency status 

%% Output:
    % betaest: (p, nv), estimators of beta(s) at grids S
    % alphaest: (p, nv), estimators of alpha(s) at grids S
    % lambdaest: (1, nv)
    % beta_ini: (p, nv), the initial values calculated by 2008 Yingcun Xia's method


[n, p] = size(x);
nv = size(y,2);
% S = rand(nv-2,1); S = sort(S); S = [0;S;1]; % Set of Sk

% tic;
%% (1) Initials by Yingcun Xia method (2008)
nd = 1;
% disp('computing initials');
% 2008 Yingcun Xia method
tmp_a = zeros(p,nv);
beta_ini = zeros(p,nv);
for s = 1:nv
    [tmp_a(:,s), ~, beta_ini(:,s), ~, ~, ~, ~, ~, ~] = dmave(x, y(:,s), nd);
end
% projection of alpha on beta
alpha_ini = zeros(p,nv);
for s = 1:nv
    %     Px_beta = beta_ini2(:,s)*(beta_ini2(:,s)'*beta_ini2(:,s))^(-1)*beta_ini2(:,s)';
    Px_beta = beta_ini(:,s)*(beta_ini(:,s))';
    alpha_ini(:,s) = tmp_a(:,s) - Px_beta * tmp_a(:,s);
end
% toc;
% smoothing the initials
betaini_smooth = zeros(p, nv);
alphaini_smooth = zeros(p, nv);
for s = 1:nv
    betaini_smooth(:,s) = locallinear0(p, S(s), 2*hs, S, beta_ini');
    alphaini_smooth(:,s) = locallinear0(p, S(s), 2*hs, S, alpha_ini');
end


%% (2) estimation at each grid by algorithm in 2016 paper
% disp('estimating coefficients at each grid seperately');
coeff_ini = [beta_ini(2:p,:); alpha_ini(2:p,:); lambda_ini];
beta2016 = zeros(p,nv);
alpha2016 = zeros(p,nv);
for s = 1:nv
    coeff = coeff_ini(:,s);
    [beta2016(:,s), alpha2016(:,s), lambdaest] = get_coeff(coeff, x, y, S, s, hx, hg, hs, 0, verbose);
end
% toc;
% smooth the estimators by algorithm in 2016 paper
beta2016_smooth = zeros(p, nv);
alpha2016_smooth = zeros(p, nv);
for s = 1:nv
    beta2016_smooth(:,s) = locallinear0(p, S(s), 2*hs, S, beta2016');
    alpha2016_smooth(:,s) = locallinear0(p, S(s), 2*hs, S, alpha2016');
end


%% (3) estimation at each grid with kernel weights
coeff_ini2 = [beta2016_smooth(2:p,:); alpha2016_smooth(2:p,:); lambda_ini];
betaest = zeros(p,nv);
alphaest = zeros(p,nv);
lambdaest = zeros(1,nv);
for s = 1:nv
    coeff = coeff_ini2(:,s); % use estimators by 2016 paper as initials
    [betaest(:,s), alphaest(:,s), lambdaest(:,s)] = get_coeff(coeff, x, y, S, s, hx, hg, hs, 1, verbose);
end
% toc;


%% (4) decide if smooth the estimators or not 
if smooth_coeff == 1
    betaest_smooth = zeros(p, nv);
    alphaest_smooth = zeros(p, nv);
    for s = 1:nv
        betaest_smooth(:,s) = locallinear0(p, S(s), 2*hs, S, betaest');
        alphaest_smooth(:,s) = locallinear0(p, S(s), 2*hs, S, alphaest');
    end
    betaest = betaest_smooth;
    alphaest = alphaest_smooth;
end


end

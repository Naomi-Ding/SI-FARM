%% Simultaneous Confidence band for beta(s) & g(x*beta(s))

function [Cb, Cg, beta_b, g_b] = SCB(x, ally, S, alphaest, betaest, lambdaest, ...
	g_est, etaest, smooth_coeff, verbose, SCB_alpha, B, idx)

% Input: 
    % smooth_coeff: 0 or 1, boolean, whether smooth_coeff the estimators (alpha(s), beta(s))
    % verbose: 0 or 1, whether to display the convergency status 
    % idx: random number generator, ranging from 1-200 for 200 replications
	% SCB_alpha: significance level, can be a vector or a scalar, the
	% length is J
	% B: number of bootstraps 

% Output: 
	% Cb: confidence bands for the coefficient functions beta(s), p x J
	% Cg: confidence band for the link function g(), scalar or a vector of
	% length J
	% beta_b: estimated beta(s) by the bootstrap samples, p x nv x B  
	% g_b: estimated g(x(beta(s)) by the bootstrap samples, n x nv x B

rng(idx);
[n, nv] = size(ally);
p = size(x,2);

std_x = std(x(:));
hx = n^(-1/3) * std_x;   %cn^(-1/3)
hg = n^(-1/5) * std_x;   %cn^(-1/5)
hs = nv^(-1/5)*0.1;

xb = x * betaest;
epsilon = ally - g_est - x * alphaest - etaest; % (n,m)


%% Bootstrap Steps
disp('Bootstrap Steps');

beta_b = zeros(p,nv,B);
g_b = zeros(n,nv,B);
dg_b = zeros(n,nv,B);

tic;
for b = 1:B
    if mod(b,10) == 0
        disp(b);
        toc;
    end
    
    %% generate random samples and construct new response
    z1 = randn(n,1);
    z2 = randn(n,nv);
    yb = g_est + x * alphaest + z1 .* etaest + z2 .* epsilon; % (n,m)
    
    %% estimate beta(s) based on new response
    %     lambda_ini = rand(1, m);
    [beta_b(:,:,b), ~, ~, ~] = coefficient_estimator(x, yb, S, hx, hg, hs, lambdaest, smooth_coeff, verbose);
    
    %% estimate g() based on new response
    Yb = yb - x * alphaest; % (n,nv)
    for i=1:n
        for s=1:nv
            xb0 = x(i, :)*betaest(:,s); %xi.T * beta(sm), scalar
            [g_b(i,s,b), dg_b(i,s,b)] = locallinear1(1, xb0, hg, xb(:), Yb(:));
        end
    end
    
end

xi_beta = repmat(betaest,1,1,B) - beta_b; % (p,m,B)
max_xi_beta = squeeze(max(abs(xi_beta), [], 2)); % (p,B)
sort_xi_beta = sort(max_xi_beta, 2); % (p,B)
Cb = sort_xi_beta(:, (1 - SCB_alpha) * B); % (p,J)

xi_g = repmat(g_est, 1,1,B) - g_b; % (n, m, B)
max_xi_g = squeeze(max(max(abs(xi_g), [],2), [], 1)); % (B,1)
sort_xi_g = sort(max_xi_g); % (B,1)
Cg = sort_xi_g( (1 - SCB_alpha) * B); % scalar or a vector of length J


end
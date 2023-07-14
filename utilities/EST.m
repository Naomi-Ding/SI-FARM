%% Estimate beta(s) & g(x*beta(s)) & the covariance function

function [betaest, alphaest, lambdaest, g_est, dg_est, etaest, beta_ini] = EST(x, ally, ...
    S, smooth_coeff, verbose, idx)
% Input: 
    % smooth_coeff: 0 or 1, boolean, whether smooth the estimators (alpha(s), beta(s))
    % verbose: 0 or 1, whether to display the convergency status 
    % idx: random number generator, ranging from 1-100 for 100 replications

% Output: 
    % betaest: estimated beta(s), p x nv 
    % alphaest: estimated alpha(s), p x nv
    % lambdaest: estimated lambda, 1 x nv 
    % g_est: estimated g(x*beta(s)), n x nv
    % dg_est: estimated first derivative of g(x*beta(s)), n x nv
    % etaest: estimated individual functions, n x nv 


    rng(idx);
    [n, nv] = size(ally);
    p = size(x,2);

    std_x = std(x(:));
    hx = n^(-1/3) * std_x;   %cn^(-1/3)
    hg = n^(-1/5) * std_x;   %cn^(-1/5)
    hs = nv^(-1/5)*0.1;

%     disp('Estimate coefficients and index function');
    %% (1) estimate coefficients beta(s) & alpha(s)
    lambda_ini = rand(1, nv);
    [betaest, alphaest, lambdaest, beta_ini] = coefficient_estimator(x, ally, S, hx, hg, hs, lambda_ini, smooth_coeff, verbose);

    %% (2) estimate g, dg
    g_est = zeros(n, nv);
    dg_est = zeros(n, nv);
    xb = x * betaest; % (n,nv)
    Y = ally - x * alphaest; % (n,nv)
    % hg = cvh1(x, betaest, xb, Y);
    for i=1:n
        for s=1:nv
            xb0 = x(i, :)*betaest(:,s); %xi.T * beta(sm), scalar
            [g_est(i,s), dg_est(i,s)] = locallinear1(1, xb0, hg, xb(:), Y(:));
        end
    end

    %% (3) Estimate covariance function Sigma_eta(s,t)
    ystar = Y - g_est; % (n, nv), ally - g(x*beta) - x*alpha
    etaest = zeros(n, nv);
    h2 = cvh2(S', ystar); % compute h2 using CV
    for i=1:n
        for s=1:nv
            etaest(i, s) = locallinear0(1, S(s), h2, S, ystar(i, :)');
        end
    end
    % Sigma_eta = etaest' * etaest / (n - 2*p); % (m,m)

end 

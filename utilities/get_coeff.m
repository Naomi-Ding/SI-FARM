%% Estimation of coefficients (beta, alpha)

function [betaest, alphaest, lambdaest] = get_coeff(coeff, x, ally, S, s, h1, h2, h, weights, verbose)

%% Input:
    % coeff: (2*p-1,1) vector, initials of (beta_{2:p}(s), alpha_{2:p}(s), lambda) for the optimization problem at grid s
    % S: (nv, 1) vector, the vector of grids for the measurement of response
    % s: scalar, ranging from 1 to nv, indicating at which grid we are estimating the coefficients
    % weights: boolean, the indicator of whether to use the kernel weights: 
            % weights == 1: use kernel weights;
            % weights ~= 1: treat each grid separately
    % verbose: whether displaying the convergency status 

%% Output:
    % betaest: (p, 1), estimators of beta(s) at grid s  
    % alphaest: (p, 1), estimators of alpha(s) at grid s
    % lambdaest: scalar


%% (1) Construct the objective function
[~,p] = size(x);
EQ_coeff = @(coeff) eq_coeff(coeff, x, ally, S, s, h1, h2, h, weights);


%% (2) Solve the optimization function using Adam algorithm 
[f, df] = EQ_coeff(coeff);
niter = 10000;
% cf0 = coeff;
% initialize
alpha = 5e-4; % step size 5e-4
beta1 = 0.9; % factor for average gradient
beta2 = 0.9999; % factor for average squared gradient
eps = 1e-8;
m = 0; v = 0; t = 1;
diff_dfnorm = -ones(niter, 1);
% while dg_norm > -3e-3 && dcf >= 1e-4 && t < niter % 9e-4
df_t = zeros(2*p-1, niter); df_t(:,t) = df;
diff_coeff = ones(niter,1);
ft = zeros(niter,1); ft(t) = f;
cf = zeros(2*p-1, niter); cf(:,t) = coeff;
while (t < niter) && (diff_dfnorm(t) < 0) ...
        && (diff_coeff(t) > 5e-4) % stoppting criteria 
    g = df_t(:,t);
    m = beta1 * m - (1 - beta1) * g;
    v = beta2 * v + (1 - beta2) * (g.^2);
    mhat = m / (1 - beta1^t);
    vhat = v / (1 - beta2^t);
    cf(:,t+1) = cf(:,t) - alpha * mhat ./ (sqrt(vhat) + eps);
    t = t + 1;
    diff_coeff(t) = norm(cf(:,t) - cf(:,t-1));
    [ft(t), df_t(:,t)] = EQ_coeff(cf(:,t));
    diff_dfnorm(t) = norm(df_t(:,t)) - norm(df_t(:,t-1));
%     cf0 = cf1;
end
% if ~isequal(verbose, 'true')
if verbose ~= 0
    fprintf('converge at %d_th iteration', t);
end


%% (3) Collect the final estimators for (beta(s), alpha(s)) at grid s
coeff_est = cf(:,t);
b_tilde = coeff_est(1:p-1); % (p-1,1)
b1 = sqrt(1 - norm(b_tilde)^2); % first component of beta(s)
betaest = [b1; b_tilde];

a_tilde = coeff_est(p:(2*p-2));
a1 = -b_tilde'*a_tilde/b1; % first component of alpha(s)
alphaest = [a1; a_tilde];

lambdaest = coeff_est(end);

end

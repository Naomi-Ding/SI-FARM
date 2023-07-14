%% Update beta(s) & alpha(s) at grid s

function [F, EQ] = eq_coeff(coeff, x, ally, S, s, hx, hg, hs, weights)

%% Input:
    % coeff: (2*p-1,1) vector, initials of (beta_{2:p}(s), alpha_{2:p}(s), lambda) for the optimization problem at grid s
    % S: (nv, 1) vector, the vector of grids for the measurement of response
    % s: scalar, ranging from 1 to nv, indicating at which grid we are estimating the coefficients
    % weights: boolean, the indicator of whether to use the kernel weights: 
            % weights == 1: use kernel weights;
            % weights ~= 1: treat each grid separately

%% Output: 
    % F: 
    % EQ: (2p-2, 1), the value of LHS of the objective function



[~,p] = size(x);
% [~,pa] = size(w);
[n,nv] = size(ally);


%% (1) Construct the intials for (beta(s), alpha(s)) at grid s 
% b_tilde = coeff(2:p); % (p-1,1)
% b1 = coeff(1);
b_tilde = coeff(1:p-1); % (p-1,1)
b1 = sqrt(1 - norm(b_tilde)^2);
beta = [b1; b_tilde];

a_tilde = coeff(p:(2*p-2));
a1 = -b_tilde'*a_tilde/b1;
alpha = [a1; a_tilde];

lambda = coeff(2*p-1);


%% (2) Compute Jacobian matrix at grid s: (beta_1,alpha_1) dependent on others
J = zeros(2*p+1,2*p-1);
J(1,1:(p-1)) = - b_tilde'/b1;
J(2:p,1:(p-1)) = eye(p-1);
J(p+1,1:(p-1)) = - a_tilde'/b1 - b_tilde'*a_tilde*b_tilde'/(b1^3);
J(p+1,p:(2*p-2)) = - b_tilde'/b1;
J((p+2):2*p,p:(2*p-2)) = eye(p-1);
J(2*p+1,2*p-1) = 1;


%% (3) Estimate g(X.T*Beta(s)) & g_dot(X.T*Beta(s)) at grid s
xb = x * beta; % (n,1)
% hg = cvh1(x, beta_ini, xb, y); % compute h1 using CV
% Y = y - w * alpha_ini; % (n,nv)
y = ally(:,s); % (n,1)
y_res = y - x * alpha; % (n,1)

g_tilde = zeros(n,1);
dg_tilde = zeros(n,1);
for i=1:n % for object i
    xb0 = x(i, :) * beta; % xi.T * beta(sm), scalar
    [g_tilde(i), dg_tilde(i)] = locallinear1(1, xb0, hg, xb(:), y_res(:));  %g(xi.T*beta(sm))
end


%% (2) Estimate alpha & beta using eq(6) at grid s
% EQ = fun_coeff(x, xb, ally, w, g_tilde, dg_tilde, alpha_ini, S, s, hx, hs); % (p,2)

%% (4) Estimate E{Xi|Xi.T*Beta(s)} & E{wi|Xi.T*Beta(s)}
mx = zeros(n,p);
% mw = zeros(n,pa);
for i = 1:n
    xb0 = xb(i); % xi.T*beta(sm), scalar
    mxi = locallinear0(p, xb0, hx, xb, x); % E[xi|xi.T*beta(sm)], (1,p)
    mx(i, :) = x(i, :) - mxi; % xi-E[xi|xi.T*beta(sm)], (1,p)
    % mwi = locallinear0(pa, xb0, hx, xb, w);  % E[wi|xi.T*beta(sm)], (1,3)
    % mw(i, :) = w(i, :) - mwi; % wi-E[wi|xi.T*beta(sm)], (1,3)
end


%% (5) Compute L.H.S of the equation, i.e., the objective function
dg_mx = dg_tilde .* mx; % (n,pb), dg(x.T*beta(s)) * {x-E[x|x.T*beta(s)]}

%% include kernel weights
if weights == 1
    all_res = ally - x * alpha - g_tilde; % (n,nv), y(sk) - w.T*alpha(s) - g(x.T*beta(s)) for all sk
    delta = S - S(s); % sk-s for all k, (nv,1)
    KH = kh(delta,hs); % Khs(sk-s) for all k, (nv,1)
    [~,nvar] = size(J);
    eq = zeros(nvar,nv);
    % NJ2s = zeros(nvar, nvar, nv);
    for sk = 1:nv
        % for sk = s
        res = all_res(:,sk); % (n,1), y(sk) - x.T*alpha(s) - g(x.T*beta(s)) at grid sk
        I1 = res .* dg_mx + 2 * lambda * repmat(beta',n,1);  % (n,p)
        I2 = res .* mx; % (n,p)
        I3 = repmat(norm(beta)^2 - 1, n,1); % (n,1)
        JI = J' * [I1, I2, I3]'; % (2p-1,n)
        JI = sum(JI,2)/n; % (2p-1,1), sum for i=1:n
        eq(:,sk) = KH(sk) * JI; % (2p-1,1)
    end
    EQ = sum(eq,2) / sum(KH); % (2p-1,1), sum for sk=1:snv
    F = sum( sum(all_res.^2,1) .* KH') / (n * sum(KH)); 
    
else
    %% treat each point seperately
    res = y_res - g_tilde; % (n,1)
    I1 = res .* dg_mx + 2 * lambda * repmat(beta',n,1);  % (n,p)
    I2 = res .* mx; % (n,p)
    I3 = repmat(norm(beta)^2 - 1, n,1); % (n,1)
    JI = J' * [I1, I2, I3]'; % (2p-1,n)
    JI = sum(JI,2)/n; % (2p-1,1), sum for i=1:n
    
    EQ = JI; % (2p,1)
    F = sum( (ally(:,s) - x*beta - g_tilde).^2) / n;
end


end
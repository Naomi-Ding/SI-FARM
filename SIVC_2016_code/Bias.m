% Bias of estimations at each grid point

% INPUT
% x : n * p matrix
% ally,ystar :  n * m matrix
% betaest : p * m matrix
% s, m, hx, hy : scalar
% t : 1 * m vevtor


% OUTPUT
% Ans : p * p matrix
% Bs : p * m matrix

% 07/25/2015
% Xinchao Luo

function [Ans, Bs] = Bias(x, ally, ystar, betaest, m, t, s, hx, hy, h)

[n, p] = size(x);
Ans = zeros(p, p);
Bs = zeros(p, n);

delta = zeros(m, 1);
for ss=1:m
    delta(ss) = t(ss) - t(s); % s - si, (m,1)
end
kernel = kh(delta, h); % weights(s-si), Kh(s-si), (m,1)
kernel = kernel/sum(kernel); % normalized, (m,1)

for ss=1:m
    xb = x*betaest(:, ss); % X.T * beta_est(s) at grid s, (n,1)
    y = ally(:, ss); % y(s) at grid s, (n,1)
    for i=1:n
        x0 = xb(i); % xi.T * beta_est(s), scalar
        [~, ddg]= locallinear1(1, x0, hy, xb, y); % g_dot(xi.T*beta_est(s))
        mx = locallinear0(p, x0, hx, xb, x); % E[xi|xi.T*beta_est(s)], (1,p)
        mx0 = x(i, :) - mx; % xi-E[xi|xi.T*beta_est(s)], (1,p)
        % weighted Seff
        Ans = Ans + kernel(ss)*ystar(i, ss)*mx0'*x(i, :)*ddg/n;
    end
end

for i=1:n
    x0 = xb(i);
    [dg, ~]= locallinear1(1, x0, hy, xb, y); % g(xi.T*beta_est(sm)) at last grid sm
    mx = locallinear0(p, x0, hx, xb, x);% E[xi|xi,T*beta_est(sm)] at last grid sm, (1,p)
    mx0 = x(i, :) - mx; % xi-E[xi|xi.T*beta_est(sm)] at last grid sm, (1,p)
    Bs(:, i) = mx0'*dg/n;
end


end
  

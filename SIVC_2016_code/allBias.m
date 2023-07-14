% Bias of all estimations

% INPUT
% x : n * p matrix
% ally,ystar :  n * m matrix
% betaest : p * m matrix
% s, m, hx, hy, h1 : scalar
% t : 1 * m vevtor


% OUTPUT
% biasbeta : p * m matrix
% biasg: n * m matrix

% 07/25/2015
% Xinchao Luo

function [biasbeta, biasg] = allBias(x, ally, ystar, betaest, m, t, hx, hy, h, h1)

[n, p] = size(x);
biasbeta = zeros(p, m);
biasg = zeros(n, m);
u2K = 0.2;

for s=1:m
    [Ans, Bs] = Bias(x, ally, ystar, betaest, m, t, s, hx, hy, h);
    tmp = zeros(p, 1);
    for i=1:n
        xb = x*betaest(:, s); % X.T * beta_est(s) at grid s, (n,1)
        y = ally(:, s); % y(s) at grid s, (n,1)
        x0 = xb(i); % xi.T * beta_est(s), scalar
        [~, ddg]= locallinear1(1, x0, hy, xb, y); % g_dot(xi.T*beta_est(s))
        biasg(i, s) = 0.5*h1^2*u2K*ddg/n;  % Thm2 (ii)
        tmp = tmp + Bs(:, i)*ddg;
    end
    biasbeta(:, s) = 0.5*h^2*Ans*tmp/n;
end

end
  

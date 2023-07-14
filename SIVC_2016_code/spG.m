% Score bootstrap for coefficients

% INPUT
% ystar, gest : n * m matrix
% x : n * p matrix
% betaest : p * m matrix
% tau : 1 * n vector
% t : 1 * m vector
% m, hx, hy, h : scalar

% OUTPUT
% sp : p * m matrix

% 11/20/2015
% Xinchao Luo

function [sp] = spG(ystar, gest,tau, x, m, t, hx, hy, h, betaest)

[p, ~] = size(betaest);
betaest_boot = zeros(p, m);
yboot = gest + ystar.*repmat(tau, 1, m);
for s=1:m
    betaest_boot(:, s) = getBeta(x, yboot, m, t, s, hx, hy, h, betaest(:, s));
end

sp = betaest_boot - betaest;

end
  

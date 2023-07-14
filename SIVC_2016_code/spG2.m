% Score bootstrap for index

% INPUT
% ystar, gest : n * m matrix
% x : n * p matrix
% betaest : p * m matrix
% tau : 1 * n vector
% m, h1 : scalar

% OUTPUT
% sp2 : n * m matrix

% 11/20/2015
% Xinchao Luo

function [sp2] = spG2(ystar, gest,tau, x, m, h1, betaest)

[n, ~] = size(x);
gest_boot = zeros(n, m);
yboot = gest + ystar.*repmat(tau, 1, m);
allxb = x*betaest;
for i=1:n
    for s=1:m
        xb0 = x(i, :)*betaest(:, s);
        gest_boot(i, s) = locallinear0(1, xb0, h1, allxb(:), yboot(:));        
    end
end

sp2 = gest_boot - gest;

end
  

% SCB
% Simultaneous confidence bands

% INPUT
% ystar, gest : n * m matrix
% x : n * p matrix
% betaest : p * m matrix
% t : 1 * m vector
% G, m, b, hx, hy, h, h1 : scalar
% betaest : p * m vector

% OUTPUT
% varCb : p * 1 vector
% varCg : scalar

% 11/10/2015
% Xinchao Luo
  
function [varCb, varCg] = SCB(ystar, gest, x, m, t, G, betaest, hx, hy, h, h1)

[n, p] = size(x);
tau = randn(n, G);
alpha = 0.05;

allCb = zeros(p, m, G);
for g=1:G
%     g
    sp = spG(ystar, gest, tau(:, g), x, m, t, hx, hy, h, betaest); 
    allCb(:, :, g) = sp;
end
varC = zeros(p, m);
for i=1:p
    for j=1:m
        tmp = allCb(i, j, :);
        tmp = sort(tmp);
        varC(i, j) = tmp(ceil(G*(1 - alpha)));
    end
end
varCb = max(varC,[], 2);

allCg = zeros(n, m, G);
for g=1:G
%     g
    sp2 = spG2(ystar, gest, tau(:, g), x, m, h1, betaest); 
    allCg(:, :, g) = sp2;
end

varC2 = zeros(n, m);
for i=1:n
    for j=1:m
        tmp = allCg(i, j, :);
        tmp = sort(tmp);
        varC2(i, j) = tmp(ceil(G*(1 - alpha)));
    end
end

varCg = max(max(varC2));

end
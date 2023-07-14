% Estiamte beta

% INPUT
% x : n * p matrix
% ally :  n * m matrix
% m, s, b, hx, hy, h : scalar
% t : 1 * m vevtor
% beta0 : p * 1 vector

% OUTPUT
% betaest : p * 1 vector

% 03/01/2015
% Xinchao Luo

function [betaest] = getBeta(x, ally, m, t, s, hx, hy, h, beta0)

% options = optimoptions('fsolve', 'Display', 'off'); 
options = optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-1, 'TolX', 1e-1); 
% fsolve(@(x) f(x), x0,option), starts at x0 to solve f(x)=0
betaest = fsolve(@(beta) sumseff(x, ally, m, t, s, beta, hx, hy, h), beta0, options);
betaest = betaest/sqrt(sum(betaest.^2));
betaest = betaest*sign(betaest(1));

end
  
% function [betaest] = getBeta(x, ally, p, m, t, b, hx, hy, beta0)
% 
% options = optimoptions('fsolve', 'Display', 'off');
% betaest = zeros(p, m);
% for s=1:m    
%     betaest(:, s) = fsolve(@(beta) meanseff(x, ally, m, t, s, beta, b, hx, hy), beta0(:, s), options);
%     betaest(:, s) = betaest(:, s)/sqrt(sum(betaest(:, s).^2));
% end
% 
% end

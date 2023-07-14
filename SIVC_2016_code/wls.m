% Sum of weighted square function 

% INPUT
% x_tra : n_tra * p matrix
% ally_tra :  n_tra * m matrix
% h2, s : scalar
% t : 1 * m vevtor
% B : p * 2 vector

% OUTPUT
% sums : scalar

% 07/26/2015
% Xinchao Luo

function [sums] = wls(x_tra, ally_tra, t, B, h2, s)

[n_tra, m] = size(ally_tra);
delta = t - t(s);
kernel = kh(delta, h2);

sums = 0;
for i=1:n_tra
    for ss = 1:m       
        xx = x_tra(i, :)*delta(ss);         
        sums = sums + (ally_tra(i, s) - x_tra(i, :)*B(:, 1) - xx*B(:, 2))^2*kernel(ss);
    end
end

end
  

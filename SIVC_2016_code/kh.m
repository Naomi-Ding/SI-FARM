% One dimensional Kernel function

% INPUT
% x : m-dimensional vector
% h : scalar

% OUTPUT
% kernel : m * 1 vector

% 02/25/2015
% Xinchao Luo

function [kernel] = kh(x, h)

tmp = x/h;
m = length(x);
kernel = zeros(m, 1);

for i=1:m
    if (abs(tmp(i))<=1)
        kernel(i) = (1 - tmp(i)^2)/h*0.75;
    end    
end
    
end

% CV for h2 

% INPUT
% t : m-dimensional vector
% ystar : n * m matrix

% OUTPUT
% h2 : scalar

% 08/25/2015
% Xinchao Luo

function [h2] = cvh2(t, ystar)

[n, m] = size(ystar);
srange = range(t);
nh = 20; CV = zeros(nh,1);
hmin = srange/m; 
hmax = srange/6; 
vh = logspace(log10(hmin),log10(hmax),nh);
etaest = zeros(n, m);

for ii=1:nh
    h = vh(ii);
    for i=1:n
        ncv = 0;
        for s=1:m
            ncv = ncv + 1;
            t_cv = t'; t_cv(ncv) = [];
            ystar_cv = ystar(i, :)'; ystar_cv(ncv) = [];
            etaest(i, s) = locallinear0(1, t(s), h, t_cv, ystar_cv);     
        end
    end
    CV(ii) = mean(mean((ystar - etaest).^2));
end

[~,flag] = min(CV);
h2 = vh(flag);

end


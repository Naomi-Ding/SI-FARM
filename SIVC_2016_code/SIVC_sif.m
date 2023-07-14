function [T_sif, h] = SIVC_sif(T, sK, ker)

if nargin < 3
    ker = @(t) exp(-0.5*t.^2);
end

[n, K] = size(T);

sK = sK(:);
s_mat = bsxfun(@minus, sK, sK');

srange = range(sK);

nh = 5; GCV = zeros(nh,1);

% crude search
hmin = srange/K; 
hmax = srange/8; 
vh=logspace(log10(hmin),log10(hmax),nh);

for ii=1:nh
    h = vh(ii);
    K_mat = feval(ker, s_mat/h)/h;
    Ks_mat = K_mat.*s_mat;
    Ks2_mat = Ks_mat.*s_mat;
    
    K_vec = sum(K_mat,2);
    Ks_vec = sum(Ks_mat,2);
    Ks2_vec = sum(Ks2_mat,2);
    
    temp_mat1 = repmat(K_vec.*Ks2_vec-Ks_vec.^2,1,K);
    temp_mat2 = repmat(Ks2_vec,1,K).*K_mat-repmat(Ks_vec,1,K).*Ks_mat;
    S_mat = temp_mat2./temp_mat1;
   
    temp1 = (eye(K)-S_mat)'*(eye(K)-S_mat);
    temp2 = (1-trace(S_mat)/K)^2;
    for nii=1:n
        GCV(ii)= GCV(ii)+T(nii,:)*temp1*T(nii,:)';
    end
    GCV(ii) = GCV(ii)/temp2;
end
[~,flag] = min(GCV);
h = vh(flag);

% delicate search
hmin = max(hmin, h-2*srange/K); 
hmax = min(hmax, h+2*srange/K); 
vh=logspace(log10(hmin),log10(hmax),nh);

for ii=1:nh
    h = vh(ii);
    K_mat = feval(ker, s_mat/h)/h;
    Ks_mat = K_mat.*s_mat;
    Ks2_mat = Ks_mat.*s_mat;
    
    K_vec = sum(K_mat,2);
    Ks_vec = sum(Ks_mat,2);
    Ks2_vec = sum(Ks2_mat,2);
    
    temp_mat1 = repmat(K_vec.*Ks2_vec-Ks_vec.^2,1,K);
    temp_mat2 = repmat(Ks2_vec,1,K).*K_mat-repmat(Ks_vec,1,K).*Ks_mat;
    S_mat = temp_mat2./temp_mat1;
   
    temp1 = (eye(K)-S_mat)'*(eye(K)-S_mat);
    temp2 = (1-trace(S_mat)/K)^2;
    for nii=1:n
        GCV(ii)= GCV(ii)+T(nii,:)*temp1*T(nii,:)';
    end
    GCV(ii) = GCV(ii)/temp2;
end
[~,flag] = min(GCV);
h = vh(flag);


K_mat = feval(ker, s_mat/h)/h;
Ks_mat = K_mat.*s_mat;
Ks2_mat = Ks_mat.*s_mat;

K_vec = sum(K_mat, 2); % K x 1 vector
Ks_vec = sum(Ks_mat, 2);
Ks2_vec = sum(Ks2_mat, 2);

temp_mat1 = repmat(K_vec.*Ks2_vec-Ks_vec.^2,1,K);
temp_mat2 = repmat(Ks2_vec,1,K).*K_mat-repmat(Ks_vec,1,K).*Ks_mat;
S_mat = temp_mat2./temp_mat1; % K x K matrix

T_sif = T*S_mat';

end






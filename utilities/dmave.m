function [gamma SEgamma B SEB G seG L seL cv] = dmave(x, y, nd);

% estimation of multiple-index model
%  Y = \gamma^\top X + G(B^\top X) + \varepsilon
%
% Input:
%     X: nxp, covariates
%     Y: nx1, response
%     nd: dimension of nonlinear part
%
% Output:
%    gamma: linear coefficient, SEgamma: SE of gamma
%    B    : dimension reduction direction in the nonlinear part.
%    SEB  : SE of B
%    G    : the link function
%    seG  : SE of G



[n,p] = size(x);
onen = ones(n,1);

mm = mean(x);
x = x - repmat(mm,n,1);
x0 = x;
ss = inv(x'*x/n + eye(p)/n/n/n*0)^0.5;
x = x*ss;


if nd == 0;
    x1 = [x ones(n,1)];
    alpha = inv(x1'*x1)*(x1'*y);
    s = std(y-x1*alpha);
    se = diag(inv(x1'*x1))*s;

    B = [];
    SEb = [];
    g = alpha(p+1);
    seG = se(p+1);
    alpha = alpha(1:p);
    SEa = se(1:p);
    cv = 0;
    for i = 1:n
        I = [1:i-1 i+1:n];
        alphai = inv(x1(I,:)'*x1(I,:))*(x1(I,:)'*y(I));
        cv = cv + (y(i) - x1(i,:)*alphai)^2;
    end
    cv = cv/n;
    L = []; seL = [];
    Chi2 = 1;
    return
end


m = p;
B = eye(p);
Ba = B;
BI = Ba;
B = B(:,1:m);
noip0 = 0;
noip1 = 1;
iterstop = 0;
Btmp = B;
rige = std(y)*mean(std(x,[],1));
%y = (y-mean(y))/std(y);
%yfit = y;

[beta B] = Mimd2(x,y,nd);  % pls simple remove this line for faster calculation
m = nd;

x0 = x;
xc = x;
enn = 1/n/n/n;

niter = floor(p*3/2);
for iter = 1:niter;
    x = xc;
    if iter >= p;
        x = x0;
    end
  	adj = p^(1/iter)*n^(1/(nd+4)-1/(m+4));
    h = p^(1/(iter+1))*mean(std(x*Ba,[],1))/n^(1/(m+4));
%    h = cvm(x*B, y-x*beta);

    h2 = 2*h*h*adj;
    
    
%    ABI = zeros(m, n*n);
    onexi = ones(n, m+1);
    for iiter = 1:4;
	dd = 0;
   	dc = 0;
        for j = 1:n;
        xij = x - repmat(x(j,:),n,1);
        sxij = sum((xij*Ba).^2,2) + sum(xij.^2,2)/1.5^iter;
        sd = sort(sxij);
      	ker = exp(-sxij/h2);  
      	ker = ker/sum(ker); 
        f = 1;
%                f = sum(exp(-sum((xij*Ba).^2,2)/h2)); 

        rker = repmat(ker, 1, p+1);
   		onexi(:,1:m) = xij*B;
   		xk = (onexi.*rker(:, 1:m+1))';
   		abi = inv(xk*onexi+eye(m+1)/n)*(xk*(y-x*beta));
      
        xdij = [x  kron(abi(1:m)', xij)];
        kxdij = (repmat(ker, 1, (m+1)*p).*xdij)';
        dd = dd + kxdij*xdij*f;
        dc = dc + kxdij*(y-abi(m+1))*f;
%        yfit(j) = abi(m+1);
        
        end

       BETA = inv(dd + eye(length(dc))*enn)*dc;
       
   	   B0 = reshape(BETA(p+1:p+p*m), p, m); %*ABI;
       [B, R] = eig(B0*B0');
       B = B(:,p-m+1:p);
       Ba = B;
       if (max(svd(B*B'-BI*BI')) < 0.001);
           break
       end
       BI = B;
       beta = (eye(p)-B*B')*BETA(1:p);       
    end

   if (max(svd(B*B'-BI*BI')) < 0.001)*(iter>p+3);
       'ok'
       break
   end
   BI = Ba;   
   
end

yfit = y;
fitcv = y;
xB = x*B;
dd = 0;
g = y;
df = y;


h2 = 2*h*h;
for j = 1:n;
      xBi = (xB - repmat(xB(j,:), n, 1));      
      sxij = sum(xBi.^2,2);
      sd = sort(sxij);
      h22 = max(h2, sd(p+1));
      ker = exp(-sxij/(h22));  
   	  onexi = [xBi/h onen];
   	  xk = (onexi.*repmat(ker, 1, m+1))';
      invxx = inv(xk*onexi+eye(m+1)*enn);
   	  abi = invxx*(xk*(y-x*beta));
      vi = invxx*(xk*x);
      fd(j,:) = abi(1:m)/h;
      vB(j,:) = vi(m+1,:)/h;
      g(j) = abi(m+1);
      df(j) = mean(ker);
      yfit(j) = abi(m+1);
%      fit0 = ker'*(y-x*beta)/(sum(ker)+1.0e-5);
       xk(:,j) = 0;
      invxx = inv(xk*onexi+eye(m+1)*enn);
   	  abi = invxx*(xk*(y-x*beta));
      fitcv(j) = abi(m+1);
%      yfit(j) = (fit0 + fit1)/2;
end
s2 = mean(df.*(y- x*beta -yfit).^2)/mean(df);
cv = mean((y- x*beta -fitcv).^2);

G = g;
seG = sqrt(0.2821*s2./(n*h^nd*df));

%subplot(2,2,2)
%plot(x*B, y-x*beta, '.')
%hold on
%plot(x*B, g, 'r.')
%subplot(2,2,1)
%plot(x*beta, y-g, '.')
%corrcoef(x*beta, x*B)

dd = x-vB;
dd0 = dd;
for i = 1:m
    dd = [dd repmat(fd(:,i),1,p).*dd0];
end
%BT = eye((m+1)*p);
%BT(1:p, 1:p) = eye(p) - B*B';
I = (df>2/(2*pi)^(nd/2)/h^nd/n);
I = df;
ddf = dd.*repmat(I, 1,size(dd,2));




dd2 = ddf'*ddf;
dd0 = dd'*ddf;


sss = kron(eye(m+1), ss);
[c d] = eig(dd0);
ddiag = diag(d);
ds = sort(ddiag);
I = find(ddiag<= ds(m*(m+1)));
ddiag(I) = 1;
ddiag = 1./ddiag;
ddiag(I) = 0;
invdd = c*diag(ddiag)*c';


S = sss*invdd*dd2*invdd*sss*s2;
SEa = sqrt(diag(S(1:p,1:p)));
alpha = ss*beta;


%%%%%%%%%%%%%%%%%%

df = y;
xa = x0*alpha;
L = y; df = y;
h = std(xa)/n^0.2;
h2 = 2*h*h;
for j = 1:n;
      xi = xa - xa(j);
      sxij = sum(xi.^2,2);
      ker = exp(-sxij/h2)/sqrt(2*pi)/h;  
      ker1 = ker.*xi;
      df(j) = mean(ker);
      c = (ker*(xi'*ker1) - ker1*sum(ker1));
      L(j) = c'*(y-g)/(sum(c)+1/n^2);
      df(j) = mean(ker);
end

seL = sqrt(0.2821*s2./(n*h*df));


%%%%%%%%%%%%%%%%%

B = ss*B;
SEb = B;
for i = 1:size(B,2);
    d = sqrt(B(:,i)'*B(:,i));
    B(:,i) = B(:,i)/d;
    A = S(i*p+1:(i+1)*p, i*p+1:(i+1)*p)/d^2;
    SEb(:,i) = sqrt(diag(A));
end
SEB = SEb;

s = inv(B'*B)^0.5;
B1 = B*s;
A = (eye(p) - 0*B1*B1');
alpha = A*alpha;
gamma = alpha;



C = [A  (-B1*kron(eye(m), alpha'))
    zeros(m*p, p) eye(m*p)];
S = C*S*C';

SEgamma = sqrt(diag(S(1:p,1:p)));


function [beta B h] = Mimd1(X,y,nd)

[n,p] = size(X);
X = X - repmat(mean(X),n,1);
S = inv(X'*X/n+eye(p)/n/n)^0.5;
X = X*S;

X1 = [X ones(n,1)];
beta = inv(X1'*X1 + eye(p+1)/n/n )*(X1'*y);
%y = y-X1*beta;

B = eye(p);
m = p;
for iter = 1:p+5;
	ab = ones(p,n);
    if iter<3
        h = n^(-1/(m+4));
    end
    if iter/3 == floor(iter/3)
        h = n^(-1/(m+4));
%        h = cvm(X*B, y);
    end
	for i = 1:n;
       	xi = X - repmat(X(i,:),n,1);
        kernel = exp(-sum((xi*B).^2,2)/(2*h*h));
        onexi = [xi ones(n,1)];
        xk = onexi.*repmat(kernel, 1, p+1);
        abi = inv(xk'*onexi+eye(p+1)*0.0001)*xk'*y;
        ab(:,i) = abi(1:p);
	end;
    ab = ab - repmat(mean(ab,2), 1, n);
	[B0 D] = eig(ab*ab');
	[D I] = sort(diag(D));
	B = B0(:,I);
    B = B(:, p+1-(1:p));
    m = max(nd, m-1);
    B = B(:,1:m);
end;

beta = S*(eye(p)-B*B')*beta(1:p);

B = S*B;
for i = 1:size(B,2);
    B(:,i) = B(:,i)/sqrt(B(:,i)'*B(:,i));
end



function [beta B h] = Mimd2(X,y,nd)

[n,p] = size(X);
X = X - repmat(mean(X),n,1);
S = inv(X'*X/n+eye(p)/n/n)^0.5;
X = X*S;

X1 = [X ones(n,1)];
beta = inv(X1'*X1 + eye(p+1)/n/n )*(X1'*y);
%y = y-X1*beta;

B = eye(p);
m = p;
for iter = 1:p+5;
	ab = ones(p,n);
    if iter<3
        h = n^(-1/(m+4));
    end
    if iter/3 == floor(iter/3)
        h = n^(-1/(m+4));
%        h = cvm(X*B, y);
    end
	for i = 1:n;
       	xi = X - repmat(X(i,:),n,1);
        kernel = exp(-sum((xi*B).^2,2)/(2*h*h));
        onexi = [xi ones(n,1)];
        xk = onexi.*repmat(kernel, 1, p+1);
        abi = inv(xk'*onexi+eye(p+1)*0.0001)*xk'*y;
        ab(:,i) = abi(1:p);
	end;
    ab = ab - repmat(mean(ab,2), 1, n);
	[B0 D] = eig(ab*ab');
	[D I] = sort(diag(D));
	B = B0(:,I);
    B = B(:, p+1-(1:p));
    m = max(nd, m-1);
    B = B(:,1:m);
end;

beta = S*(eye(p)-B*B')*beta(1:p);

B = S*B;
for i = 1:size(B,2);
    B(:,i) = B(:,i)/sqrt(B(:,i)'*B(:,i));
end



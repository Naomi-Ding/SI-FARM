%% simu_data_gen: the function for generating the simulation datasets
% Input: 
    % n: sample size,
    % nv: number of grids 
    % nsimu: number of replications
    % G0: a q x p matrix, the linear relationship between x & z
    % ul: severity of heterogeneity, an integer chosen from [0,1,2,3] 
    % index_g: choice of the index function g(x*beta(s))
    % seed: the seed for the random number generator 
% Output:
    % beta0: ground truth of beta(s) measured at the grids S
    % alpha0: ground truth of alpha(s)
    % alpha_star: ground truth of G' * alpha(s)
    % S: a vector of length nv, the grids
    % simu_x: generated x, n x p x nsimu
    % simu_y: generated response, n x nv x nsimu 
    % simu_eta: generated individual functions, n x nv x nsimu 
    % simu_eta_star: simu_eta + epsilon*alpha_star, n x nv x nsimu
    % simu_g: simulated index function, n x nv x nsimu
    % simu_z: simulated unobserved confounders, n x q x nsimu


function [beta0, alpha0, alpha_star, S, simu_x, simu_y, simu_eta, ...
    simu_eta_star, simu_g, simu_z] = simu_data_gen(n, nv, nismu, ... 
        G0, ul, index_g, sigma_xi, sigma_epsilon, seed)

    rng(seed);

    % Scenarios on ul
    if ul == 0        % (i) ul = 0, no heterogeneity
        u = 0;
    elseif ul == 1    % (ii) ul ~ U(0,0.2), weak heterogeneity
        u = 0.2 * rand;
    elseif ul == 2    % (iii) ul ~ U(0.3,0.5), moderate heterogeneity
        u = 0.3 + 0.2 * rand;
    elseif ul == 3    % (iv) ul ~ U(0.6,0.9), high heterogeneity
        u = 0.6 + 0.3 * rand;
    end
    G = u .* G0;
    [q,p] = size(G);

    % sigma_xi = [0.5 0.5]/5;
    % sigma_epsilon = 0.05;
    S = rand(nv-2,1); S = sort(S); S = [0;S;1]; % Set of Sk

    % beta(s) & alpha(s)
    beta0 = zeros(p, nv);
    for s = 1:nv
        sk = S(s);
        beta0(:,s) = [1 + sk^2; 4 * sk * (1-sk); (1- sk^2 + 2*sk)];
        beta0(:,s) = beta0(:,s) / sqrt(sum(beta0(:,s).^2));
    end
    alpha0 = (S'.^3 + 1)/4;
    alpha_star = G' * alpha0; % (p,nv)


    %% Generate Data
    % g & psi function
    psi1 = @(s) 2*s - 1; % eigen functions for eta
    psi2 = @(s) 1;
    if index_g == 1
        g = @(xb) sin(2*xb) + 2 * cos(2 + xb);
        dg = @(xb) 2 * cos(2*xb) - 2 * sin(2 + xb);
    elseif index_g == 2
        g = @(xb) exp(-xb);
        dg = @(xb) - exp(-xb);
    end

    % Simulate x
    simu_x = randn(n,p,nsimu);
    simu_x = (simu_x - mean(simu_x,1)) ./ std(simu_x,1);

    % Simulate z & eta & epsilon
    simu_eta = zeros(n,nv,nsimu);
    simu_eta_star = zeros(n,nv,nsimu);
    simu_g = zeros(n,nv,nsimu);
    simu_y = zeros(n,nv,nsimu);
    simu_z = zeros(n,q,nsimu); % zi = wi.T*gamma + ei, ei~N(0,0.1^2)
    simu_e = 0.01*randn(n,q,nsimu);
    for nn = 1:nsimu    
        simu_z(:,:,nn) = simu_x(:,:,nn) * G' + simu_e(:,:,nn); % 
        simu_eta(:,:,nn) = sigma_xi(1)*randn(n,1) * psi1(S)' + ...
            sigma_xi(2)*randn(n,1) * psi2(S)'; %(n,nv)
        simu_eta_star(:,:,nn) = simu_e(:,:,nn)*alpha0 + simu_eta(:,:,nn); % (n,nv)
        for s = 1:nv
            xb = simu_x(:,:,nn) * beta0(:,s); % X.T*beta(sk), (n,1)
            simu_g(:,s,nn) = g(xb); % (n,1)
            xa = simu_x(:,:,nn) * alpha_star(:,s); % X.T*alpha*(sk), (n,1)
            simu_y(:,s,nn) = simu_g(:,s,nn) + xa + simu_eta_star(:,s,nn) + ...
                sigma_epsilon*randn(n,1);
        end
    end


end 
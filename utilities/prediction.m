%% prediction: calculate the prediction errors (squared error & absolute error) 
            % for PLSIM & SIVC & MVCM respectively for one dataset
            % the dataset is split to a training set and a test set, 
            % with the proportion pp, generally pp=0.3, 0.5 or 0.7 

% Input:
    % smooth_coeff: 0 or 1, boolean, whether smooth the estimators (alpha(s), beta(s))
    % verbose: 0 or 1, whether to display the convergency status 
	% pp: proportion of training data & test data
	% idx: seed for random number generator, ranging from 1-100 for 100 replications
% Output: 
	% SE: a vector of length 3, containing the squared prediction error 
	       % for PLSIM, SIVC & MVCM respectively
	% AE: a vector of length 3, containing the absolute prediction error 
	       % for PLSIM, SIVC & MVCM respectively


function [SE, AE] = prediction(x, ally, S, smooth_coeff, verbose, pp, idx)

	rng(idx);

	n_tra = floor(n*pp);
	n_tes = n - n_tra;
    ind = 1:n;
    ind_tra = randsample(n, n_tra);
    ind_tes = ind;
    ind_tes(ind_tra) = [];
    
    x_tra = x(ind_tra, :); ally_tra = ally(ind_tra, :);
    x_tes = x(ind_tes, :); ally_tes = ally(ind_tes, :);

    SE = zeros(3, 1);
	AE = zeros(3, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLSIM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Prediction for PLSIM')
    % Estimation    
    [betaest, alphaest, lambdaest, g_est, dg_est, etaest] = EST(x_tra, ally_tra,...
        S, smooth_coeff, verbose, nn);
    
    % Prediction
    % disp('computing response using test data');
    tes_g = zeros(n_tes, nv);
    allxb = x_tes * betaest;
    Ytes = ally_tes - x_tes * alphaest;
    for i=1:n_tes
        for s=1:nv
            xb0 = x_tes(i, :)*betaest(:, s);
            tes_g(i, s) = locallinear0(1, xb0, hg, allxb(:), Ytes(:));
        end
    end
    ally_tes_est = tes_g + x_tes * alphaest;
    SE(1) = (sum(sum((ally_tes_est - ally_tes).^2))/n_tes/nv)^0.5;
    AE(1) = (sum(sum(abs(ally_tes_est - ally_tes)))/n_tes/nv);    
    toc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SIVC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Prediction for SIVC');
    % Computing initials by 2008 Yingcun Xia method
    nd = 1;
    beta_ini = zeros(p,nv);
    for s = 1:nv
        [~, ~, beta_ini(:,s), ~, ~, ~, ~, ~, ~] = dmave(x, ally(:,s), nd);
    end

    betaest_SIVC = zeros(p, nv);
    for s=1:nv
        betaest_SIVC(:, s) = getBeta(x_tra, ally_tra, nv, S', s, hx, hg, h, beta_ini(:, s));
    end
    
    ally_tes_est_SIVC = zeros(n_tes, nv);
    allxb = x_tes*betaest_SIVC;
    for i=1:n_tes
        for s=1:nv
            xb0 = x_tes(i, :)*betaest_SIVC(:, s);
            ally_tes_est_SIVC(i, s) = locallinear0(1, xb0, hg, allxb(:), ally_tes(:));
        end
    end
    SE(2) = (sum(sum((ally_tes_est_SIVC - ally_tes).^2))/n_tes/nv)^0.5;
    AE(2) = (sum(sum(abs(ally_tes_est_SIVC - ally_tes)))/n_tes/nv);
    toc;
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MVCM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('prediction for MVCM');
    ally_tes_est_MVCM = zeros(n_tes, nv);
    betaest_MVCM = zeros(p, nv);
    for s=1:nv
        options = optimoptions('fminunc', 'Algorithm','Quasi-Newton', 'Display', 'off');
        B0 = zeros(p, 2);
        B0(:, 1) = betaest_SIVC(:, s); B0(:, 2) = 0;
        Best = fminunc(@(B) wls(x_tra, ally_tra, S', B, h2, s), B0, options);
        betaest_MVCM(:, s) = Best(:, 1);
    end

    for s=1:nv
        ally_tes_est_MVCM(:, s) = x_tes* betaest_MVCM(:, s);
    end
    SE(3) = (sum(sum((ally_tes_est_MVCM - ally_tes).^2))/n_tes/nv)^0.5;
    AE(3) = (sum(sum(abs(ally_tes_est_MVCM - ally_tes)))/n_tes/nv);
    toc;

end 

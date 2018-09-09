% Estimates the means (mus, lambdas, and gammas) and the variances (sigmas and tau) given W.

function [mus_hat,sigmas_hat,deltas_hat,gammas_hat,tau_hat] = TCA_fit_model_means_vars(X,W,mus_hat,sigmas_hat,tau_hat,C2,deltas_hat,C1,gammas_hat)

%% Paramters for the optimization
MAX_ITERATIONS = 5;
MIN_VAR = 10^-7;
EPSILON = 0.001; % percent of difference for convergence of optimization

n = size(X,1);
m = size(X,2);
k = size(W,2);
p1 = size(C1,2);
p2 = size(C2,2);

% Settings for the optimization of mus, gammas and lambdas
options_means = optimset('lsqlin');
options_means = optimset(options_means,'Display','none');
lb = ones(1,k+p2+k*p1)*-inf;
ub = ones(1,k+p2+k*p1)*inf;
lb(1:k) = 0;
ub(1:k) = 1;
% Settings for the optimization of the sigma_{jk} values
options_vars = optimset('fmincon');
options_vars = optimset(options_vars,'Display','none','GradObj','on');
lb_vars_k = ones(k,1)*MIN_VAR;
lb_vars = ones(k+1,1)*MIN_VAR;

C1_= create_interactions_matrix(W,C1);

% If the parameters are all zeros, need to find initial estimates first.
if (sum(sum((mus_hat == 0))) == m*k)
    
    % Get an initial estimate of mus, deltas, gammas under the assumptions
    % tau=0 and sigmas_{1j},...,sigmas_{kj} for each j.
    W_norms = sqrt(sum(W.^2,2));
    X_tilde = X./repmat(W_norms,1,m);
    W_tilde = W./repmat(W_norms,1,k);
    C1_tilde = C1_./repmat(W_norms,1,p1*k);
    C2_tilde = C2./repmat(W_norms,1,p2);
    for j = 1:m
        x = lsqlin([double(W_tilde) C2_tilde C1_tilde],double(X_tilde(:,j)),[],[],[],[],lb,ub,[],options_means);
        mus_hat(j,:) = x(1:k)';
        deltas_hat(j,1:p2) = x(k+1:k+p2)';
        gammas_hat(j,1:k*p1) = x(k+p2+1:k+p2+k*p1)';
    end
    
    sigmas_hat = repmat(std(X)'./k,1,k);   % use for a starting point for the optimization
    
end

% Iteratively improve the estimates of the means and the variances.
global_ll_new = -minus_log_likelihood_model_tau(X,W,mus_hat,sigmas_hat,tau_hat,C2,deltas_hat,C1_,gammas_hat);
for iter = 1:MAX_ITERATIONS
    
    global_ll_current = global_ll_new;
      
    % for each sigma_j, estimate sigma_{1j},...,sigma_{kj},tau (i.e. initially assume
    % that tau is specific for each site j).
    for j = 1:m
        x0 = zeros(k+1,1);
        x0(1:k) = sigmas_hat(j,:)';
        x0(k+1) = tau_hat;
        f = @(x) minus_log_likelihood_site(X(:,j),W,mus_hat(j,:)',x(1:k),x(k+1),C2,deltas_hat(j,:)',C1_,gammas_hat(j,:)',1);
        [x,~,~,~,~,~] = fmincon(f,double(x0),[],[],[],[],lb_vars,[],[],options_vars);
        sigmas_hat(j,:) = x(1:k)';
    end
    
    % Estimate tau.
    f = @(x) minus_log_likelihood_model_tau(X,W,mus_hat,sigmas_hat,x,C2,deltas_hat,C1_,gammas_hat);
    [x,~,~,~,~,~] = fmincon(f,double(tau_hat),[],[],[],[],MIN_VAR,[],[],options_vars);
    tau_hat = x;
    
    % for each sigma_j, estimate sigma_{1j},...,sigma_{kj}
    for j = 1:m
        x0 = sigmas_hat(j,:)';
        f = @(x) minus_log_likelihood_site(X(:,j),W,mus_hat(j,:)',x,tau_hat,C2,deltas_hat(j,:)',C1_,gammas_hat(j,:)',0);
        [x,~,~,~,~,~] = fmincon(f,double(x0),[],[],[],[],lb_vars_k,[],[],options_vars);
        sigmas_hat(j,:) = x';
    end
            
    % Estimate mus, deltas, gammas
    for j = 1:m
        W_norms = sqrt( sum((W.^2).*repmat(sigmas_hat(j,:).^2,n,1),2)+tau_hat^2 );
        X_j_tilde = X(:,j)./W_norms;
        W_tilde = W./repmat(W_norms,1,k);
        C1_tilde = C1_./repmat(W_norms,1,p1*k);
        C2_tilde = C2./repmat(W_norms,1,p2);
        x = lsqlin([double(W_tilde) C2_tilde C1_tilde],double(X_j_tilde),[],[],[],[],lb,ub,[],options_means);
        mus_hat(j,:) = x(1:k)';
        deltas_hat(j,1:p2) = x(k+1:k+p2)';
        gammas_hat(j,1:k*p1) = x(k+p2+1:k+p2+k*p1)';
    end
    
    global_ll_new = -minus_log_likelihood_model_tau(X,W,mus_hat,sigmas_hat,tau_hat,C2,deltas_hat,C1_,gammas_hat);
    global_ll_diff = global_ll_new-global_ll_current;
    if (global_ll_diff < EPSILON*abs(global_ll_current))
        break
    end

end

end
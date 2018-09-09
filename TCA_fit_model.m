% Fits the TCA model.
% 
% INPUT:
% X - an n by m data matrix, where n is the number of individuals and m is
% the number of sites.
% W = an n by k matrix of cell proportion estimates for each
% individual, where k is the number of assumed cell types.
% C1 = an n by p1 matrix of p1 cell-type specific covariates for each
% individual. If there are no cell-type specific covariates pass [].
% C2 = an n by p2 matrix of p2 global covariates for each individual. If there are no global covariates pass [].
% fit_W (optional) = 1 / 0 for refiting W (default is 0)

function [W,mus_hat,sigmas_hat,tau_hat,deltas_hat,gammas_hat] = TCA_fit_model(X, W, C1, C2, fit_W)

% Parameters for the optimization
MAX_ITERATIONS = 10;
EPSILON = 0.001; % percent of difference for convergence of optimization;

if (~exist('fit_W','var'))
    fit_W = 0;
end
assert(fit_W == 0 | fit_W == 1);

disp('Starting TCA_fit_model...')

n = size(X,1);
m = size(X,2);
k = size(W,2);
p1 = size(C1,2);
p2 = size(C2,2);

% init parameters
mus_hat = zeros(m,k);
sigmas_hat = zeros(m,k);
deltas_hat = zeros(m,p2);
gammas_hat = zeros(m,p1*k);
tau_hat = 0.001;

if (p1 == 0)
    C1 = zeros(n,0);
end
if (p2 == 0)
    C2 = zeros(n,0);
end

disp('iteration 0');
% Learn initial estimates of the means and variances.
[mus_hat,sigmas_hat,deltas_hat,gammas_hat,tau_hat] = TCA_fit_model_means_vars(X,W,mus_hat,sigmas_hat,tau_hat,C2,deltas_hat,C1,gammas_hat);

% Iteratively improve the estimates of W and the estimated of the means and the variances.
C1_= create_interactions_matrix(W,C1);
global_ll_new = -minus_log_likelihood_model_tau(X,W,mus_hat,sigmas_hat,tau_hat,C2,deltas_hat,C1_,gammas_hat);

if (fit_W)
    for iter = 1:MAX_ITERATIONS

        disp(['iteration ' num2str(iter)])
        global_ll_current = global_ll_new;

        W = TCA_fit_model_W(X,W,mus_hat,sigmas_hat,tau_hat,C2,deltas_hat,C1,gammas_hat);
        C1_= create_interactions_matrix(W,C1);
        [mus_hat,sigmas_hat,deltas_hat,gammas_hat,tau_hat] = TCA_fit_model_means_vars(X,W,mus_hat,sigmas_hat,tau_hat,C2,deltas_hat,C1,gammas_hat);

        global_ll_new = -minus_log_likelihood_model_tau(X,W,mus_hat,sigmas_hat,tau_hat,C2,deltas_hat,C1_,gammas_hat);
        global_ll_diff = global_ll_new-global_ll_current;
        disp(['Improvement in global log likelihood: ' num2str(global_ll_diff)])
        if (global_ll_diff < EPSILON*abs(global_ll_current))
            break
        end

    end
end

disp('TCA_fit_model is done.')

end

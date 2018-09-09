% Fits the model of the conditional distribution Y|X_j.
%
% INPUT:
% y - an n-length phenotype vector (n is the number
% of individuals).
% W - an n by k matrix of cell proportions for each
% individual, where k is the number of assumed cell types.
% X_j - an n-length vector of methylation levels in site j.
% mus_j, sigmas_j - k length vectors of the cell-type-specific means and
% standard deviations in site j.
% tau - the standard deviation of the i.i.d. component of variation in the
% model.
% C1 = an n by p1 matrix of p1 cell-type specific covariates for each
% individual. If there are no cell-type specific covariates pass [].
% C2 = an n by p2 matrix of p2 global covariates for each individual. If there are no global covariates pass [].
% C3 = an n by p3 matrix of p3 covariates to account for in the association
% test with the phenotype (these may be the same as the phenotypes in C1
% and C2). If there are no covariates pass [].
% deltas_j - a p2-length vector of the effect sizes of the global covariates C2.
% gammas_j - a p1*k length vector of the effect sizes of the
% cell-type-specific covariates C2.
% cell_types - vector with the cell types (i.e. their indices) to be tested for association with y.
% betas_0 - vector of length length(cell_types) with initial guess of the
% effect sizes of the cell types in cell_types on the phenotype
% alphas_0 - p3-length vector with initial guess of the
% effect sizes of the covariates in C3.
% phi_0 - initial guess of the i.i.d. component of variation in the conditional model.
% single_beta - 0/1 (default is 0). If 1 then the optimization constrain
% all the cell types in cell_types to have the same effect size beta.

function [phi,beta,alpha] = fit_conditional_model(y, W, X_j, mus_j, sigmas_j, tau, C1, C2, C3, deltas_j, gammas_j, cell_types, betas_0, alphas_0, phi_0, single_beta)

if (~exist('single_beta','var'))
    single_beta = 0;
end

VAR_MIN = 0.00001;
l = length(cell_types);

% initial value for the optimization.
if (single_beta)
    x0 = [phi_0 betas_0(1) alphas_0]';
else
    x0 = [phi_0 betas_0 alphas_0]';
end
% define the function for minimization
if (single_beta)
    fun = @(x) conditional_model_minus_log_likelihood([x(1); repmat(x(2),length(cell_types),1); x(3:end)], y, W, X_j, mus_j, sigmas_j, tau, C1, C2, C3, deltas_j, gammas_j, cell_types);
else
    fun = @(x) conditional_model_minus_log_likelihood(x, y, W, X_j, mus_j, sigmas_j, tau, C1, C2, C3, deltas_j, gammas_j, cell_types);
end

lb = -Inf*ones(length(x0),1);
lb(1) = VAR_MIN;
ub = Inf*ones(length(x0),1);

if (single_beta)
    options = optimoptions(@fmincon,'Display','none','GradObj','off');
    l = 1;
else
    options = optimoptions(@fmincon,'Display','none','GradObj','on');
end

x = fmincon(fun,x0,[],[],[],[],lb,ub,[],options);
phi = x(1);
beta = x(2:(l+1));
alpha = x((l+2):end);

if (single_beta)
    beta = repmat(beta,length(cell_types),1);
end

end

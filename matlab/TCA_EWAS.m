% Performs EWAS under the TCA model with cell-type-specific effects.
%
% INPUT:
% y - an n by 1 phenotype vector
% X - an n times m matrix of observed methylation levels for n individuals
% at m sites
% W - an n times k matrix of cell-type proportions
% mus_hat - an m times k matrix of cell-type specific means (given by
% TCA_fit_model).
% sigmas_hat - an m times k matrix of cell-type specific standard
% deviations (given by TCA_fit_model).
% tau_hat - estimated standard deviation of the measurement noise (given by
% TCA_fit_model).
% C1 (optional) - an n times p1 matrix of p1 cell-type-specific covariates
% (default is []); make sure to use the same C1 used in TCA_fit_model.
% C2 (optional) - an n times p2 matrix of p2 global covariates (default is
% []); make sure to use the same C1 used in TCA_fit_model.
% deltas_hat (optional) - an m times p2 matrix of effect sizes for the p2 global covariates for
% each of the sites (default is []); (given by TCA_fit_model)
% gammas_est (optional) - an m times p1*k matrix of cell-type-specific effect sizes for the p1 covariates for
% each of the sites  (default is []); (given by TCA_fit_model)
% C3 (optional) - an n times p3 matrix of p3 covariates that may affect the phenotype; do not include an intercept term (default is []).
% cell_types (optional) - a vector with the cell types to test for cell-type specific
% methylation associations with the phenotype y (default is all cell
% types in a joint test, i.e. 1:k).
% single_beta (optional) - an indicator (i.e. 0 or 1). If 1 then the
% optimization constrains all the cell types in cell_types to have the same
% effect size beta (default is 0).
%
% OUTPUT:
% pvals - pvalues (for all m sites) for the effect sizes of the cell-type
% specific methylation levels of the cell types in cell_types. 
% betas_hat - an m times length(cell_types) matrix of effect sizes for the
% cell types in cell_types (for all m sites).
% alphas_hat - an m times (p3+1) matrix with the estimated effect sizes of the
% p3 covariates in C3 (the first value is an intercept).
% methylation levels of the cell types in cell_types.
% phis_hat - an m-length vector with the estimates of phi for each site.
% ll1s - an m-length vector with the log likelihood of the (alternative)
% TCA model for each site.

function [pvals,betas_hat,alphas_hat,phis_hat,ll1s] = TCA_EWAS(y, X, W, mus_hat, sigmas_hat, tau_hat, C1, C2, deltas_hat, gammas_hat, C3, cell_types, single_beta)

n = length(y);
m = size(X,2);
k = size(W,2);

if (~exist('single_beta','var'))
    single_beta = 0;
end

if (~exist('C1','var') | isempty(C1))
    C1 = zeros(n,0);
    gammas_hat = zeros(m,0);
end
if (~exist('C2','var') | isempty(C2))
    C2 = zeros(n,0);
    deltas_hat = zeros(m,0);
end

if (~exist('C3','var') | isempty(C3))
    C3 = zeros(n,0);
end
C3 = [ones(n,1) C3];
p3 = size(C3,2);

if (~exist('cell_types','var'))
    cell_types = 1:k;
end

% Likelihood of the null model
mdl0 = fitlm(C3(:,2:end),y);
ll0 = mdl0.LogLikelihood;
phi_0 = sqrt(mdl0.SSE/n);
alphas_0 = mdl0.Coefficients.Estimate';
betas_0 = zeros(1,length(cell_types));

pvals = zeros(m,1);
ll1s = zeros(m,1);
betas_hat = zeros(m,length(cell_types));
alphas_hat = zeros(m,p3);
phis_hat = zeros(m,1);
for j = 1:m
    if (mod(j,1000) == 0)
        disp(['site ' num2str(j) ' out of ' num2str(m) ' sites...']);
        disp(['current median and min: ' num2str(median(pvals(1:j))) ', ' num2str(min(pvals(1:j))) ]);
    end
    % Try using the estimates of phi,alpha given by the null model as the
    % starting point for the optimization.
    [phi1,beta1,alpha1] = fit_conditional_model(y, W, X(:,j), mus_hat(j,:)', sigmas_hat(j,:)', tau_hat, C1, C2, C3, deltas_hat(j,:)', gammas_hat(j,:)', cell_types, betas_0, alphas_0, phi_0, single_beta);
    x = [phi1; beta1; alpha1];
    ll1 = -conditional_model_minus_log_likelihood(x, y, W, X(:,j), mus_hat(j,:)', sigmas_hat(j,:)', tau_hat, C1, C2, C3, deltas_hat(j,:)', gammas_hat(j,:)', cell_types);
    ll1s(j) = ll1;
    if (single_beta)
        pvals(j) = chi2cdf(2*(ll1-ll0),1,'upper');
    else
        pvals(j) = chi2cdf(2*(ll1-ll0),length(cell_types),'upper');
    end
    betas_hat(j,:) = beta1';
    alphas_hat(j,:) = alpha1';
    phis_hat(j) = phi1;
end

end
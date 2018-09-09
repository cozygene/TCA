% Calculates the minus log likelihood of the model and the corresponding gradient.
% 
% INPUT:
% x - (1+length(cell_types)+p3)-length column vector, which includes the following parameters: phi, betas (length(cell_types) values), and alphas (p3 effect sizes for the p3 covariates in C3).
% y - an n-length column vector of phenotypic levels
% W - an n times k matrix of cell type proportions
% X_j - an n-length column vector of methylation levels at site j
% mus_j - a k-length column vector of cell-type specific means for site j
% sigmas_j - a k-length column vector of cell-type specific standard deviations for site j
% tau - the standard deviation of the i.i.d. component in the model
% C1 = an n by p1 matrix of p1 cell-type specific covariates for each
% individual. If there are no cell-type specific covariates pass
% zeros(n,0);
% C2 = an n by p2 matrix of p2 global covariates for each individual. If there are no global covariates pass zeros(n,0).
% C3 = an n by p3 matrix of p3 covariates to account for in the association
% test with the phenotype (these may be the same as the phenotypes in C1
% and C2). If there are no such covariates pass zeros(n,0).
% deltas_j - a p2-length vector of the effect sizes of the global
% covariates C2 in site j.
% gammas_j - a p1*k length vector of the effect sizes of the
% cell-type-specific covariates C2 in site j.
% cell_types - vector with the cell types (i.e. their indices) to be tested for association with y.

function [ll,g] = conditional_model_minus_log_likelihood(x, y, W, X_j, mus_j, sigmas_j, tau, C1, C2, C3, deltas_j, gammas_j, cell_types)

k = length(mus_j);
l = length(cell_types);
phi = x(1);
betas = x(2:l+1);
alphas = x(l+2:end);
n = length(y);
p1 = size(C1,2);

X_j_tilde = X_j - C2*deltas_j - W*mus_j - sum(W.*(C1*reshape(gammas_j,p1,k)),2);
t0 = tau^2 + W.^2*sigmas_j.^2;
t1 = W(:,cell_types)*((sigmas_j(cell_types).^2).*betas);
t2 = t1./t0;
t3 = t1.*t2;
t4 = phi^2 + (betas.^2)'*sigmas_j(cell_types).^2 - t3;
gamma_indices = zeros(p1*l,1);
for i = 1:l
	gamma_indices(1+(i-1)*p1:p1*i) = 1+(cell_types(i)-1)*p1:cell_types(i)*p1;
end
t5 = repmat(mus_j(cell_types),1,n)' + C1*reshape(gammas_j(gamma_indices),p1,l) + repmat(X_j_tilde,1,l).*W(:,cell_types).*(repmat(sigmas_j(cell_types).^2,1,n)')./repmat(t0,1,l);
t6 = C3*alphas + t5*betas;
t7 = t6-y;
t8 = t7.^2;

% Calculate the negative log likelihood
ll = 0.5*( n*log(2*pi) + sum(log(t4)) + sum(t8./t4));

% Calculate the gradient of the negative log likelihood
g = zeros(length(x),1);
g(1) = phi*( sum(1./t4) - sum(t8./(t4.^2)) );
g(2:l+1) = (sigmas_j(cell_types).^2)'.*sum( (repmat(betas,1,n)' - W(:,cell_types).*repmat(t2,1,l))./repmat(t4,1,l) ) + sum( (repmat(mus_j(cell_types),1,n)' + C1*reshape(gammas_j(gamma_indices),p1,l) + W(:,cell_types).*(repmat(sigmas_j(cell_types).^2,1,n)').*repmat(X_j_tilde./t0,1,l)).*repmat(t7./t4,1,l) ) -(sigmas_j(cell_types).^2)'.*sum( (repmat(betas,1,n)'-W(:,cell_types).*repmat(t2,1,l)).*repmat(t8./(t4.^2),1,l) );
g(l+2:end) = [C3'*(t7./t4)]';


end
% The negative of the log likelihood of the TCA model and its gradient in
% site j.
% INPUT:
% X_j - a column vector with the levels across samples in site j
% W - n times k cell composition matrix
% mus_j a column vector of the cell-type-specifi means in site j
% tau
% C2 - n times p2 global covariates (pass zeros(n,0) if no covariates)
% deltas_j - a column vector of p2 coefficients for the p2 global
% covariates
% C1 - n times p1 cell-type-specific covariates (pass zeros(n,0) if no covariates)
% gammas_j - a column vector of p1*k coefficients for the p1
% cell-type-specific covariates (in each cell type)
% eval_tau_grad - an indicator variable (0/1), indicating whether to
% evaluate an (k-1)-length gradient (where the (k+1)-th component is the
% negative of the derivative with respect to tau).
function [fval,grad] = minus_log_likelihood_site(X_j,W,mus_j,sigmas_j,tau,C2,deltas_j,C1_,gammas_j,eval_tau_grad)
    
n = length(X_j);
k = size(W,2);
V = (W.^2)*(sigmas_j.^2)+tau^2;
U = (W*mus_j + C2*deltas_j + C1_*gammas_j - X_j).^2;
fval = -double( -0.5*(n*log(2*pi) + sum(log(V))) - 0.5*sum(U./V) );
grad = -double( sum( ((W.^2).*repmat(sigmas_j',n,1).*repmat(U,1,k))./repmat(V.^2,1,k) ) - sum( ((W.^2).*repmat(sigmas_j',n,1))./repmat(V,1,k) ) )';

if (eval_tau_grad)
    g = zeros(k+1,1);
    g(k+1) = -(tau*(sum( U./(V.^2)) - sum(1./V)));
    g(1:k) = grad;
    grad = g;
end

end
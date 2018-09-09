% The negative of the log likelihood of the TCA model and its derivative with respect to tau.
% INPUT:
% X - n times m matrix
% W - n times k cell composition matrix
% mus - m times k matrix of cell-type-specific means
% tau
% C2 - n times p2 global covariates (pass zeros(n,0) if no covariates)
% deltas - m times p2 matrix of coefficients for the p2 global covariates
% C1 - n times p1 cell-type-specific covariates (pass zeros(n,0) if no covariates)
% gammas_j - m times p1*k matrix of coefficients for the p1
% cell-type-specific covariates (in each cell type)
function [fval,grad] = minus_log_likelihood_model_tau(X,W,mus,sigmas,tau,C2,deltas,C1_,gammas)
    
n = size(X,1);
m = size(X,2);

fval = 0;
grad = 0;
for j = 1:m
    sigmas_j = sigmas(j,:)';
    mus_j = mus(j,:)';
    deltas_j = deltas(j,:)';
    gammas_j = gammas(j,:)';
    X_j = X(:,j);
    V = sum((W.^2).*repmat(sigmas_j.^2,1,n)',2)+tau^2;
    U = (sum(W.*repmat(mus_j,1,n)',2) + C2*deltas_j + C1_*gammas_j - X_j).^2;
    f = -double( -0.5*(n*log(2*pi) + sum(log(V))) - 0.5*sum(U./V) );
    g = -double( sum(tau*U./(V.^2)) - tau*sum(1./V));
    fval = fval + f;
    grad = grad + g;
end

end
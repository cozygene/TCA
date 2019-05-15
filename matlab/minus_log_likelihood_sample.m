% The negative of the log likelihood of the TCA model and its gradient in
% sample i.
% INPUT:
% X_j - a column vector with the levels in smaple i across all site
% W_i - a k length column vector of cell composition of sample i
% mus - m times m matrix of cell-type-specific means
% sigmas - m times m matrix of cell-type-specific sigmas
% tau
% C2 - p2 length column vector of global covariates of sample i (pass zeros(1,0) if no covariates)
% deltas - m times p2 matrix of coefficients for the p2 global covariates
% C1_i_ - p1 length column vector of cell-type-specific covariates for sample i (pass zeros(1,0) if no covariates)
% gammas - m times p1*k matrix of coefficients for the p1 cell-type-specific covariates (in each cell type)
function [fval,grad] = minus_log_likelihood_sample(X_i,W_i,mus,sigmas,tau,C2_i,deltas,C1_i,gammas)
    
m = length(X_i);
k = length(W_i);
p1 = length(C1_i);

C1_i_ = create_interactions_matrix(W_i',C1_i')';

V = sum((sigmas.^2).*repmat(W_i.^2,1,m)',2)+tau^2;
U = deltas*C2_i + gammas*C1_i_ + sum(mus.*repmat(W_i,1,m)',2) - X_i;

fval = -double( -0.5*(m*(log(pi*2)) + sum(log(V)) + sum((U.^2)./V)) );

C_tilde = zeros(m,k);
for h = 1:k
    C_tilde(:,h) = gammas(:,1+(h-1)*p1:h*p1)*C1_i;
end
grad = -double( -(sum((repmat(W_i,1,m)'.*sigmas.^2)./repmat(V,1,k)) + sum(( (mus + C_tilde).*repmat(U,1,k).*repmat(V,1,k) - (repmat(W_i,1,m)'.*sigmas.^2).*repmat(U.^2,1,k) )./(repmat(V.^2,1,k)))) )';

end
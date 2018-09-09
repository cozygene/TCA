% Gets the paramters of the Tensor Composition Analysis model and estimates
% the Z_{hj}^i values.
%
% INPUT:
% X - n by m data matrix (n samples, m features)
% W - n by k weights matrix
% mus - m by k matrix of mean levels for each feature in each source
% sigmas- m by k matrix of variances for each feature in each source
% C2 - n by p2 matrix of "global" covariates (i.e. covariates that do not have source-specific effect size)
% deltas - m by p1 matrix of effect sizes for the global covariates
% C1 - n by p1 matrix of  source-specific covariates
% gammas - m by k*p2 matrix of effect sizes for the source-specific
% covariates; each p1 consecutive columns are the effect sizes of the p2
% covariates in one particular source (the first p2 items for the first
% source, the next p2 iterms for the second soruce and so on).
%
% OUTPUT:
% Z - cell array of size k with an n by m matrix for each of the k sources.

function Z = TCA(X, W, mus, sigmas, tau, C2, deltas, C1, gammas)

n = size(X,1);
m = size(X,2);
k = size(W,2);
C3 = [C1 C2];
Z = cell(k,1);
p1 = size(C1,2);
p2 = size(C2,2);
for h = 1:k
    Z{h} = zeros(n,m);
end
U = 0;
for j = 1:m
    Sig = diag(sigmas(j,:).^2);
    Sig_inv = inv(Sig);
    for i = 1:n
        if ((p1+p2) > 0)
            U = Sig_inv * [repmat(deltas(j,:)',1,k) ; reshape(gammas(j,:),p1,k)]' * C3(i,:)';
        end
        a_ij = inv((W(i,:)'*W(i,:))./(tau^2) + Sig_inv) * ((X(i,j)/(tau^2))*W(i,:)' + Sig_inv*mus(j,:)' + U);
        for h = 1:k
            Z{h}(i,j) = a_ij(h);
        end
    end
end

end

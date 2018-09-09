% Gets W, an n by k matrix of cell composition and C1, an n times p1 matrix
% of cell-type specific covariates
% Creates C1_, a matrix with interaction terms of the cell-type specific
% covaraites and the cell-type composition.
function C1_= create_interactions_matrix(W,C1)

n = size(W,1);
k = size(W,2);
p1 = size(C1,2);
W_rep = zeros(n,p1*k);
for h = 1:k
	W_rep(:,1+p1*(h-1):p1*h) = repmat(W(:,h),1,p1);
end
C1_ = repmat(C1,1,k) .* W_rep;

end

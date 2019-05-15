% Fits the conditional model, assuming the effect size beta is 0 (i.e. the null model).
% y - an n length vector of phenotype
% C - an n times p matrix of covaraites, potentially affecting y. C should contain at least one column of constant values.
% Returns phi, the estimated standard deviation, and alpha, the estimated effect sizes of the covariates.

function [phi, alpha] = fit_conditional_null_model(y, C)
	alpha = regress(y, C);
	phi = sqrt(sum(((C*alpha - y).^2))/(length(y)-1));
end

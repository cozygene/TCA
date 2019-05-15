% Estimates the means (mus, lambdas, and gammas) and the variances (sigmas and tau) given W.

function W = TCA_fit_model_W(X,W,mus_hat,sigmas_hat,tau_hat,C2,deltas_hat,C1,gammas_hat)

n = size(X,1);
k = size(W,2);

% Settings for the optimization of W
options_W = optimset('fmincon');
options_W = optimset(options_W,'Display','none','GradObj','on');
lb_W = zeros(k,1);
ub_W = ones(k,1);

% Estimate w_{1i},...,w_{ki} for each smaple i
Aeq = ones(1,k);
for i = 1:n
    f = @(x) minus_log_likelihood_sample(X(i,:)',x,mus_hat,sigmas_hat,tau_hat,C2(i,:)',deltas_hat,C1(i,:)',gammas_hat);
    [x,~,~,~,~,~] = fmincon(f,double(W(i,:)'),[],[],Aeq,1,lb_W,ub_W,[],options_W);
    W(i,:) = x';
end

end
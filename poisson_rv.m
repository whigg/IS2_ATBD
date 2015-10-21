function x=poisson_rv(lambda, N)

% want mult*lambda >>50 -> mult=50/lambda
% want 1/mult < 1 -> mult > 1
mult=max(10,50/lambda);
x=round(sum(rand(ceil(mult*lambda),N)<(1/mult)));

 
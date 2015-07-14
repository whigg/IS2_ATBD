function x=poisson_rv(lambda, N)

% simple rejection-rule generator for Poisson random variables.
% Generates a set of M * lambda uniform random numbers (0<r<1), 
% counts the number that are < 1/M.  M is chosen so that M lambda >>
% lambda for large lambda, and so that the round-off error in the 
% calculation is adequately small for small lambda.

if lambda==0; 
    x=zeros(1,N);
    return
end

% want mult*lambda >>50 -> mult=50/lambda
% want 1/mult < 1 -> mult > 1
mult=max(10,min(100, 50/lambda));
x=sum(rand(ceil(mult*lambda),N)<(1/mult));

 
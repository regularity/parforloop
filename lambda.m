% -------------------------------------------------------------------------
%  lambda.m
%
%  Compute the optimal hazard function, L(q), given Chebyschev coefficients 
%  for the value function V(q), vc, the maximizing value q*, and the mean
%  frequency of reviews, l
% -------------------------------------------------------------------------



function lq = lambda(qval,l,vc,qstar,qspace,theta,kappa)

[nq,cq] = size(qval); 

V = @(x) funeval(vc,qspace,x);

Max = V(qstar).*ones(nq,cq);
K = kappa.*ones(nq,cq);
Val = zeros(nq,cq); 

parfor cols=1:cq
    Val(:,cols)= V(qval(:,cols));
end

if theta==0
    lq = double((Max-Val-K)>0); 
    
else
    if theta==inf
        exponent = zeros(size(Max));
        
    else
        
        exponent = (1/theta).*(Max-Val-K);
        
    end
    
    lq1 = (l/(1-l)).*exp(exponent);
    
    intermed = ones(size(exponent)) + lq1;
    
    invintermed = intermed.^(-1);
    
    lq = ones(size(exponent))-invintermed; 
    
end


return


 

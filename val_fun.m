% -------------------------------------------------------------------------
%  val_fun.m
%
%  Given the fixed point problem V(q) = T(V(q)), and a set of Chebyshev 
%  coefficients, vc, for the value function V(q), this function 
%  computes the residual between V(q) and T(V(q)), evaluated at the set  
%  of collocation nodes, (qval). 
% 
%  To be used in Broyden's method of function approximation, which finds 
%  Chebyshev coefficients that minimize this residual up to a specified
%  tolerance level. 
%
% Called by Get_Policy.m (with vbroyden)
% Calls maxeval (with newton)
% -------------------------------------------------------------------------


function [resid,qstar,auxVal] = ...
         val_fun(vc,vbas,vbas_new,qval,qstar,lbar,s,fs,beta,kappa,theta,PI,cutoffs)

     
% Initializations:

nq = size(qval,1); 
ns = size(s,1); 
qmin = min(qval);
qmax = max(qval);

  
% Intermediate expressions: 

% current value function (based on current guess for vc coeffs):
vfn = vbas*vc;                      % nq:1

% max value function (based on current guess for vc coeffs):
qstar = newton('max_eval',qstar,vc,qmin,qmax);
bas_star = chebbas(nq,qmin,qmax,qstar); 
vfnstar = bas_star*vc;

% value if update (based on current guess for vc coeffs):
E0 = vfnstar-kappa;

% next state value function (based on current guess for vc coeffs):
vfnew = vbas_new*vc;
EE  = ones(size(vfnew)); 


% Algorithm varies by theta:
if theta==cutoffs(1,1)                          % theta = 0
    
    E = max(E0.*EE,vfnew);
    
elseif theta>cutoffs(1,1) && theta<cutoffs(2,1) % low theta (baseline case)
    
    E1 = theta.*log((1-lbar).*exp((1/theta).*(vfnew-E0.*EE))+lbar.*EE);
    E1mat = reshape(E1,nq,ns);
    
    E0mat = E0.*ones(size(E1mat));
    
    E = E0mat+E1mat; 
  
elseif theta>=cutoffs(2,1) && theta<cutoffs(3,1) % high theta
    
    E = (1-lbar).*vfnew+lbar.*E0.*EE+(1/2).*(1/theta).*lbar.*(1-lbar).*((E0.*EE-vfnew).^2);
    
elseif theta==cutoffs(3,1)                       % theta = infinity
    
    E = (1-lbar).*vfnew+lbar.*E0.*EE;
    
end   


% New value function:
Emat = reshape(E,nq,ns);  
S = Emat*fs;
auxVal = PI+beta.*S;


% Update return vector: 
resid = vfn-auxVal;    
    
return  



 

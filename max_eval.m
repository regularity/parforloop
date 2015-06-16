% -------------------------------------------------------------------------
%  max_eval.m
%
%  Given a set of Chebyschev coefficients, vc, and a current estimate of
%  the maximizing value, q0, this function evaluates the criterion 
%  function and its Jacobian. 
%
%  To be used in Newton's method for computing q*(vc).
% -------------------------------------------------------------------------

function [mval,mjac] = max_eval(q0,vc,qmin,qmax)

nv=size(vc,1);

% Compute criterion function: 
bas1=chebbas(nv,qmin,qmax,q0,1);
mval=bas1*vc;

% Compute its derivative: 
bas2=chebbas(nv,qmin,qmax,q0,2);
mjac=bas2*vc;

return 


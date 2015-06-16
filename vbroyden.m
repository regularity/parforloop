% VBROYDEN Computes root of function via Broyden's Inverse Method
%
% Adapted from Toolbox broyden.m to also return maximizing value
% Uses valfun.m
%
% Uses inverse Jacobian estimate and backstepping
% USAGE
%   [x,fval,fjacinv] = broyden(f,x,varargin)
% INPUTS
%   f       : name of function of form:
%               fval=f(x,optional additional parameters)
%   x       : initial guess for root (d by 1)
%   varargin: additional arguments for f [optional]
% OUTPUTS
%   x       : root of f (d by 1)
%   fval    : function value estimate (d by 1)
%   fjacinv : inverse Jacobian estimate (d by d)
%   qstar   : value of q that maximizes value function for given x 
%
% Setable options (use OPTSET):
%   maxit     : maximum number of iterations
%   tol       : convergence tolerance
%   maxsteps  : maximum number of backsteps
%   showiters : display results of each iteration
%   initb     : an initial inverse Jacobian aprroximation matrix
%   initi     : if initb is empty, use the identity matrix to initialize
%               if 0, a numerical Jacobian will be used

% Copyright (c) 1997-2002, Paul L. Fackler & Mario J. Miranda
% paul_fackler@ncsu.edu, miranda.4@osu.edu

function [x,qstar,fval,fjacinv] = vbroyden(f,x,varargin)

maxit     = optget('broyden','maxit',100);
tol       = optget('broyden','tol',sqrt(eps));
maxsteps  = optget('broyden','maxsteps',25);
showiters = optget('broyden','showiters',0);
initb     = optget('broyden','initb',[]);
initi     = optget('broyden','initi',0);

if maxsteps<1, maxsteps=1; end
if isempty(initb)
  if initi
    fjacinv=eye(size(x,1));
  else
    fjacinv = fdjac(f,x,varargin{:});
    fjacinv = inv(fjacinv);  
  end
else
  fjacinv = initb;
end

[fval,qstar] = feval(f,x,varargin{:});
fnorm = norm(fval,'inf');
for it=1:maxit
   if fnorm<tol, return; end
   dx = -(fjacinv*fval);
   fnormold = inf;
   for backstep=1:maxsteps
      [fvalnew,qstar] = feval(f,x+dx,varargin{:});
      fnormnew = norm(fvalnew,'inf');
      if fnormnew<fnorm, break, end
      if fnormold<fnormnew
        fvalnew=fvalold;
        fnormnew=fnormold; 
        dx=dx*2; 
        break
      end
      fvalold  = fvalnew;
      fnormold = fnormnew;
      dx = dx/2;
   end
   x = x+dx;
   if any(isnan(x)|isinf(x))
     error('Infinities or NaNs encountered.')
   end

   if fnormnew>fnorm
     if initi
       fjacinv=eye(size(x,1));
     else 
       fjacinv = fdjac(f,x,varargin{:});
       fjacinv = inv(fjacinv);
     end
   else
      temp = fjacinv*(fvalnew-fval);
      fjacinv = fjacinv + (dx-temp)*(dx'*fjacinv/(dx'*temp));
   end
   fval=fvalnew;
   fnorm=fnormnew;
   if showiters, fprintf('%4i %4i %6.2e\n',[it backstep fnorm]); end
end
warning('Failure to converge in broyden');
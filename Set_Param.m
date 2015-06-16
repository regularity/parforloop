% -------------------------------------------------------------------------
% Set_Param.m
% Program for Labor. Sets paramters and settings for algorithm. 
% Last update: June 2015. 
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Parameters 
% -------------------------------------------------------------------------

periods     = 12;                               % Frequency of model (12 = monthly)
rate_annual = 0.03;                             % Steady state real interest rate
delta       = (1/(1+rate_annual))^(1/periods);  % Discount factor
epsi        = 5;                                % Elasticity of substitution
eta         = 1/epsi;                           % Inverse of elasticity of substitution
alpha       = 1;                                % Exponent on employment in production function 
beta        = 1.06;                             % Exponent on hours in production function 
nu          = 0.9;                              % Inverse of Frisch elasticity
gamma       = ((1-eta)*(alpha*(1+nu)-beta))/((1+nu)-beta*(1-eta)) ;         % Curvature of profit function 
gamma2      = (1-eta)/(1-alpha*(1-eta));                                    % Scaler of innovation in employment gap law of motion
if gamma <=0 || gamma >=1, error('The parameter gamma must be between 0 and 1.'), end

abar        = 0;                                % Mean of (inverse of) productivity shocks
asig        = (10).*1.e-2/sqrt(periods);      % Standard deviation of productivity shocks <--------------- TO BE REVIEWED
avar        = asig^2;                           % Variance of productivity shocks
arho        = 1;                                % Serial correlation of productivity  
na          = 75;                               % Discretization of productivity shock
[a,ha]      = qnwnorm(na,abar,avar);            % ha approximates pdf of a~N(abar,avar)

kappa       = 0.1;                              % Fixed cost of undertaking a review
theta       = 0.6;                              % Unit cost of review signal


% -------------------------------------------------------------------------
% Settings
% -------------------------------------------------------------------------

fprintf(' Settings .... ');

% Profit function 
pf   = @(z) exp(gamma.*z) - gamma.*exp(z) ; 
% dpf  = @(z) gamma.*exp(gamma.*z) - gamma.*exp(z) ; 
% ddpf = @(z) (gamma^2).*exp(gamma.*z) - gamma.*exp(z) ; 


% Innovations to employment gap 
s           = gamma2.*a;                        % Gridpoints
sbar        = gamma2*abar;                      % Mean of innovation to gap
ssig        = gamma2*asig;                      % Standard deviation of innovation to gap
svar        = ssig^2;                           % Variance of innovation to gap
ns          = length(s);        
hs          = ha;      


% Discretization of gap z
nz          = 101;                              % Number of grid points 
zdispersion = 50;                               % Dispersion for innovation to gap
zmax        = zdispersion*ssig;                 % Range for gaps 
zmax        = round(zmax*100)./100 ;            % DO I NEED THIS ? 
zmin        = -zmax; 


% Optimal Cheb nodes used in Cheb approximations of V(z), J1(z), JL(z)
vspace = fundefn('cheb',nz,zmin,zmax);          % Cheb polynomial of (nz-1)degree  
z      = funnode(vspace);                       % Cheb nodes 
zmin   = min(z); 
zmax   = max(z); 


% Next state (for value function approximation) 
s1                  = kron(ones(size(s)),z);
s2                  = kron(s,ones(size(z)));
znew                = s1+s2;   
znew((znew<zmin),1) = zmin;
znew((znew>zmax),1) = zmax;  
hs2                 = kron(hs,ones(size(z)));   


% Algorithm settings 
tol_lbar     = 1e-8;                                 % Tolerance for lbar 
max_it_lbar  = 30;                                   % Maximum number of iterations for lbar
cutoffs      = [0;15;inf];                           % Algorithm differs across theta of 0, low, high, infinity: 
plz          = 1;                                    % Print results every lbar iterations (0 for no display)
options      = optimset('Display','off','TolX', 1.0e-10, 'TolFun',1.0e-10);


% Simulation settings
T            = 500;                                  % Number of periods to compute everything
seed         = 10;                                   % Seed for simulations in PriceSimulations.m
l_series     = 10*12;                                % Length of each simulated series 10 years each year has 12 months
n_series     = 1000;                                 % Number of simulated series per set of parameters. 
baseline     = 1;                                    % Baseline parameterization <--------------- TO BE REVIEWED


% Initializations
PI           = pf(z);                                % Per-period profit 
Voz          = PI;                                   % Initial guess: value fn
vc           = funfitxy(vspace,z,Voz);               % Initial guess: Cheb coeffs
vbas         = chebbas(nz,zmin,zmax,z);              % Vfn = vbas*vc
vbas_new     = chebbas(nz,zmin,zmax,znew);           % Vfn_new = vbas_new*vc

zstar        = 0;                                    % Initial guess: reset value upon review 
high         = 0.9; low = 0; 
guess        = @(x,y) 1/2*(x+y); 
lbar         = guess(low,high);                      % Initial guess: lambda-bar

Jlambda      = zeros(max_it_lbar,1);
Jstar        = zeros(max_it_lbar,1);

J1z0         = (delta/(1-delta*(1-lbar))).*ones(size(z));      % Initial guess: J0
JLz0         = (delta*lbar/(1-delta*(1-lbar))).*ones(size(z)); % Initial guess: J1
jc1          = funfitxy(vspace,z,J1z0);
jcL          = funfitxy(vspace,z,JLz0);


% Functions for future updating
V   = @(x,c) funeval(c,vspace,x);
L   = @(x,l,c,zs) lambda(x,l,c,zs,vspace,theta,kappa);     
J1  = @(x,j1) funeval(j1,vspace,x);
JL  = @(x,jL) funeval(jL,vspace,x);








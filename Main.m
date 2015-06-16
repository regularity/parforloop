% -------------------------------------------------------------------------
% MAIN.m
%
% Main program for Labor
% Calls Set_Param, Get_Policy, Run_Sim
% Requires Miranda and Fackler's CompEcon Toolbox      % <------------------------------- ADD TO PATH !
% Code based on ICSDP
%
% Last update: June 2015 
% -------------------------------------------------------------------------

clear
close all
clc
format compact

fid      = fopen('Output.txt','w');
filename = ['mats/Results_', datestr(date), '.mat'] ; 

% Add COMPECON toolbox
cepath='C:\My Computer\G-Dropbox\Dropbox\DunRA\RA6_parcode\labor1_ic_par\compecon\';
path([cepath 'CEtools'],path); path([cepath 'CEdemos'],path)

% Setup Parfor Pool that connects to 4 nodes
myPool = parpool(4);

% -------------------------------------------------------------------------
% Run baseline
% -------------------------------------------------------------------------

RUN = 1; 
Set_Param; 
Get_Policy; 
 
if exist(filename,'file')==2
    load(filename)
    
else
    theta_g{RUN}     = theta;
    kappa_g{RUN}     = kappa;
    gamma_g{RUN}     = gamma;
    asig_g{RUN}      = asig; 
    ssig_g{RUN}      = ssig; 
    lbar_g{RUN}      = lbar;
    zstar_g{RUN}     = zstar;
    z_g{RUN}         = z; 
    vc_g{RUN}        = vc; 
    V_g{RUN}         = V; 
    L_g{RUN}         = L; 
end

% -------------------------------------------------------------------------
% Changing theta
% -------------------------------------------------------------------------

% Unpack selected kappa 
kappa = kappa_g{baseline};

% Range of thetas
fprintf(fid, ['\n Policy iterations for kappa = ',num2str(kappa),' and theta = [0.16 - 160].', '.... \n']);
type('Output.txt')
range_theta = [0.16 1.6 16 160] ;

% Policy iterations 
Set_Param;
for theta = range_theta
    
    RUN = RUN + 1; 
    
    % Get policy
    Get_Policy; 
    
    % Save solution 
    theta_g{RUN}     = theta;
    kappa_g{RUN}     = kappa;
    gamma_g{RUN}     = gamma;
    asig_g{RUN}      = asig;  
    ssig_g{RUN}      = ssig; 
    lbar_g{RUN}      = lbar;
    zstar_g{RUN}     = zstar;
    z_g{RUN}         = z; 
    vc_g{RUN}        = vc; 
    V_g{RUN}         = V; 
    L_g{RUN}         = L; 
    
    % Store solution 
    save(filename, '*_g'); 
      
    % Reset 
    Set_Param;
    % close all
end

% Compare hazard functions: 

% set(0,'defaultaxescolororder',[0 0 0; 0.5 0.5 0.5]) % black and gray
set(0,'DefaultAxesColorOrder',[0 0 0]) % only black 
set(0,'DefaultAxesLineStyleOrder',{'-','--','-.',':','-o','-*'}) 

fig_hazard_thetas = figure; 
hold all 

parfor ii = 1:length(range_theta)
    
    plot(z, L_g{ii+1}(z,lbar_g{ii+1}, vc_g{ii+1}, zstar_g{ii+1})) ; % ,plot_style{ii}) ;  % ii+1 because RUN = 1 is the baseline parameterization 
    leg_text_thetas{ii} = ['$\theta = $',num2str(theta_g{ii+1})]; 
    
end

xlabel('$z$','Interpreter','latex')
ylabel('$\Lambda(z)$','Interpreter','latex')

leg_fig = legend(leg_text_thetas, 'Location', 'Southeast'); 
set(leg_fig,'Interpreter','latex','FontSize',12);
legend boxoff
 
saveas(fig_hazard_thetas, 'figures/fig_hazard_thetas', 'pdf'); 
saveas(fig_hazard_thetas, 'figures/fig_hazard_thetas', 'eps'); 


% -------------------------------------------------------------------------
% Changing kappa
% -------------------------------------------------------------------------

% Unpack selected theta
theta = theta_g{baseline};

% Range of kappas
fprintf(fid,['\n Policy iterations for theta = ',num2str(theta),' and kappa = [0.05 - 0.45].', '.... \n']);
type('Output.txt')
range_kappa = [0.05 0.15 0.25 0.35 0.45] ;

% Policy iterations 
Set_Param;
for kappa = range_kappa
    
    RUN = RUN + 1; 
    
    % Get policy
    Get_Policy; 
    
    % Save solution 
    theta_g{RUN}     = theta;
    kappa_g{RUN}     = kappa;
    gamma_g{RUN}     = gamma;
    asig_g{RUN}      = asig;  
    ssig_g{RUN}      = ssig; 
    lbar_g{RUN}      = lbar;
    zstar_g{RUN}     = zstar;
    z_g{RUN}         = z; 
    vc_g{RUN}        = vc; 
    V_g{RUN}         = V; 
    L_g{RUN}         = L; 
    
    % Store solution 
    save(filename, '*_g'); 
      
    % Reset 
    Set_Param;
    % close all

end

% Compare hazard functions 

fig_hazard_kappas = figure; 
hold all 

parfor ii = 1:length(range_kappa)
    
    plot(z, L_g{end-ii+1}(z,lbar_g{end-ii+1}, vc_g{end-ii+1}, zstar_g{end-ii+1})) ;  % ii+1 because RUN = 1 is the baseline parameterization 
    leg_text_kappas{ii} = ['$\kappa = $',num2str(kappa_g{end-ii+1})]; 
    
end

xlabel('$z$','Interpreter','latex')
ylabel('$\Lambda(z)$','Interpreter','latex')

leg_fig = legend(leg_text_kappas, 'Location', 'Southeast'); 
set(leg_fig,'Interpreter','latex','FontSize',12);
legend boxoff
 
saveas(fig_hazard_kappas, 'figures/fig_hazard_kappas', 'pdf'); 
saveas(fig_hazard_kappas, 'figures/fig_hazard_kappas', 'eps'); 


fclose('all');

% % To remove color and line style settings (without restarting Matlab): 
% set(0, 'DefaultAxesColorOrder', 'remove')
% set(0,'DefaultAxesLineStyleOrder','remove') 


return 


% -------------------------------------------------------------------------
% Run simulations 
% -------------------------------------------------------------------------

Run_Sim;    
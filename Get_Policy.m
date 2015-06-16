% -------------------------------------------------------------------------
%  Get_Policy.m
%
% Program for Labor. Iterates to convergence on lambda-bar
% Calls valfun, solve_fp. Uses vbroyden, broyden, and other MF Toolbox
% functions
% Last update: June 2015 
% -------------------------------------------------------------------------


%fprintf(fid,'\n Finding lambda_bar .... ');
%type('Output.txt')

fprintf('\n Finding lambda_bar .... ');

ldiff  = 0;
lit   = 1;

while lit < max_it_lbar
    tic
    lit = lit+1; 
    
    % Get value function and reset value corresponding to conjectured lbar
    [vc,zstar]  = vbroyden('val_fun',vc,vbas,vbas_new,z,zstar,lbar,s,hs,delta,kappa,theta,PI,cutoffs);   %   Q: why is val_fun using lbar and not lambda? (5/16/2012)
    Jstar(lit) = zstar; 

    % Get updated hazard function
    Lz_new = L(znew, lbar, vc, zstar);
    
    % Obtain new estimate of lbar
    Lzmat = reshape(Lz_new,nz,ns);
    D1    = beta;
    DL    = beta.*(Lzmat*hs);
    jc1   = broyden('solve_fp',jc1,vbas,vbas_new,z,Lz_new,D1,s,hs,delta);
    jcL   = broyden('solve_fp',jcL,vbas,vbas_new,z,Lz_new,DL,s,hs,delta);
    Jlambda(lit) = JL(zstar,jcL)/J1(zstar,jc1);
    
    % Check convergence
    ldiff = Jlambda(lit)-lbar;
    if abs(ldiff) < tol_lbar, break, end;
    
    % Update guess
    if ldiff > 0
        low  = lbar;
    elseif ldiff < 0
        high = lbar;
    end
    
    lbar = guess(low,high);
    
    % Print status
    if plz>0 && mod(lit,plz)==0
        clc
        fprintf('\n ******************************************************')
        fprintf(['\n * Finding lambda-bar in RUN ', num2str(RUN)])        
        fprintf(['\n *       lbar diff: ',num2str(ldiff)])
        fprintf(['\n *       lbar iteration: ',num2str(lit)])
        fprintf('\n ******************************************************\n')
    end
    toc       
    
end


% Display output 
fprintf(fid,' ******************************************************')      %#ok<*PRTCAL>
fprintf(fid,['\n Solution for RUN ' num2str(RUN), ' with kappa = ',num2str(kappa),', theta = ',num2str(theta),' found']);
fprintf(fid,'\n * Output:')
fprintf(fid,['\n *       lbar: ',num2str(lbar)])
fprintf(fid,['\n *       zstar: ',num2str(zstar)])
fprintf(fid,'\n * Convergence')        
fprintf(fid,['\n *       lbar diff: ',num2str(ldiff)])
fprintf(fid,['\n *       lbar iteration: ',num2str(lit)])
fprintf(fid,'\n * Input:')
fprintf(fid,['\n *       theta: ',num2str(theta)])
fprintf(fid,['\n *       kappa: ',num2str(kappa)])
fprintf(fid,['\n *       gamma: ',num2str(gamma)])
fprintf(fid,['\n *       asig: ',num2str(asig)])
fprintf(fid,['\n *       ssig: ',num2str(ssig)])
fprintf(fid,'\n * Settings:')
fprintf(fid,['\n *       na: ',num2str(na)])
fprintf(fid,['\n *       nz: ',num2str(nz)])
fprintf(fid,['\n *       zdispersion: ',num2str(zdispersion)])
fprintf(fid,'\n ******************************************************\n')
type('Output.txt')

% Plot profit function
fig_profit=figure; plot(z,pf(z), 'b-')
    xlabel('$z$','Interpreter','latex')
    leg_text = '$\pi(z)$'; 
    leg_fig = legend(leg_text, 'Location', 'Southwest'); 
    set(leg_fig,'Interpreter','latex','FontSize',12);
    legend boxoff
    
  %  xlim([-1.5 1.5])
    saveas(fig_profit,['figures/fig_profit_',num2str(RUN)],'pdf')   
    saveas(fig_profit,['figures/fig_profit_',num2str(RUN)],'eps')  

% Plot hazard function
fig_hazard=figure; plot(z,L(z,lbar,vc,zstar), 'b-')
    xlabel('$z$','Interpreter','latex')
    leg_text = '$\Lambda(z)$'; 
    leg_fig = legend(leg_text, 'Location', 'Northwest'); 
    set(leg_fig,'Interpreter','latex','FontSize',12);
    legend boxoff
    
   % xlim([-1.5 1.5])
    saveas(fig_hazard,['figures/fig_hazard_',num2str(RUN)],'pdf')   
    saveas(fig_hazard,['figures/fig_hazard_',num2str(RUN)],'eps')  






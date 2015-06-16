% -------------------------------------------------------------------------
%  solve_fp.m
%
%  Solves contraction mapping of the form 
%      F(q) = C(q) + beta*E{[1-L(q_next)]*F(q_next)}
% 
%  Called by Get_Policy (with broyden) to solve for J10(q), JL0(q) 
%            and in getLinL to solve for Vl(q), J1l(q) and JLl(q).
%
%  Last update: October 2010, January 2011. 
% -------------------------------------------------------------------------


function [resid,auxFl]=solve_fp(cl,vbas,vbas_new,qval,Lq_new,Cq,s,fs,beta)

nq = size(qval,1); 
ns = size(s,1); 

Fl = vbas*cl;                       % f_current (nq:1)
Fl_new = vbas_new*cl;               % f_next (nq*na:1)

E = (1-Lq_new).*Fl_new; 
Emat = reshape(E,nq,ns);

S = Emat*fs; 

auxFl = Cq + beta.*S;

resid = Fl-auxFl;                   % return vector to be min'd

return 
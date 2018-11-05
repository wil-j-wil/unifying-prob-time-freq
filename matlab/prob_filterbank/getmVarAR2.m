function mVar = getmVarAR2(Lam,Var)

% function mVar = getmVarAR2(Lam,Var)
%
% Computes the marginal variance of an AR2 process
% i.e. <x_t^2>
%
% INPUTS
% Lam = dynamical parameters, size [D,2] (D = number of AR2 processes)
% Var = innovations variance, size [D,1]
%
% OUTPUTS
% mVar = marginal variances, size [D,1] 

lam1 = Lam(:,1);
lam2 = Lam(:,2);

mVar = (1-lam2).*Var./(1-lam1.^2-lam2.^2-lam2-lam1.^2.*lam2+lam2.^3);
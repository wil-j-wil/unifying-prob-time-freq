function [Lam,Var] = cosFreq2AR2(cosCF,cosDF,mVar);
  
% function [Lam,Var] = cosFreq2AR2(cosCF,cosDF,mVar);
%
% Converts cosine centre frequency and cosine bandwidth into AR2
% parameters.
%
% INPUTS
% cosCF = cos centre frequencies size [D,1]
% cosDF = cos bandwidths size [D,1]
%
% OUTPUTS
% Lam = dynamical parameters, size [D,2] (D = number of AR2 processes)
% Var = innovations variance, size [D,1]

a1 = cosDF.^2+4*cosCF.^2;
a2 = -2*cosDF.^2+8*cosCF.^2-2;
z2 = -a1/2-real(sqrt(a1.^2/4-a2+2));
Lam(:,2) = (z2+real(sqrt(z2.^2-4)))/2;
Lam(:,1) = 4*cosCF.*Lam(:,2)./(Lam(:,2)-1);
Var = mVar.*(1-Lam(:,1).^2-Lam(:,2).^2-Lam(:,2)-Lam(:,1).^2.* ...
	     Lam(:,2)+Lam(:,2).^3)./(1-Lam(:,2));
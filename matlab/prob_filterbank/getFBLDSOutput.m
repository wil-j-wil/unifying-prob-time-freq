function [S,covS] = getFBLDSOutput(Xfin,Pfin)

% function [S,covS] = getFBLDSOutput(Xfin,Pfin);
%
% Get the output in the correct format from the FB Kalman Smoother.
%
% INPUTS
% Xfin = mean of the latent variables, size [N,2D,T] (N=1)
% Pfin = covariance of the latent variables, size [2D,2D,T]
%
% OUTPUTS
% S = mean of the latent variables, size [D,T]
% covS = covariance of the latent variables, size [2D,2D,T]

[N,TwoD,T] = size(Xfin);  

D =  TwoD/2;

indRe = [1:2:TwoD-1];
indIm = [2:2:TwoD];

S = squeeze(Xfin(1,indRe,:)+i*Xfin(1,indIm,:));

covS = zeros(TwoD,TwoD,T);

covS = [Pfin(indRe,indRe,:),Pfin(indRe,indIm,:); 
	Pfin(indIm,indRe,:),Pfin(indIm,indIm,:)];


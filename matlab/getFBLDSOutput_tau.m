function [S,varargout] = getFBLDSOutput_tau(Xfin,Pfin,tau)

% function [S,covS] = getFBLDSOutput(Xfin,Pfin);
%
% Get the output in the correct format from the FB Kalman Smoother.
%
% INPUTS
% Xfin = mean of the latent variables, size [N,2D,T] (N=1)
% Pfin = covariance of the latent variables, size [2D,2D,T]
%
% OUTPUTS
% S = mean of the latent variables (just the observable dimensions), size [D,T]
% OPTIONAL OUTPUTS
% covS = covariance of the latent variables (just the observable dimensions), size [2D,2D,T]
% Sfull = mean of the latent variables (everything), size [tauD,T]
% covSfull = covariance of the latent variables (everything), size [2tauD,2tauD,T]

varargout = cell(1,nargout-1);

[N,TwoDtau,T] = size(Xfin);

if nargout > 2
    
    indRe = [1:2:TwoDtau-1];
    indIm = [2:2:TwoDtau];

    Sfull = squeeze(Xfin(1,indRe,:)+i*Xfin(1,indIm,:));
    
    S = Sfull(1:tau:end,:);
    if nargout > 3
        covSfull = [Pfin(indRe,indRe,:),Pfin(indRe,indIm,:); 
                    Pfin(indIm,indRe,:),Pfin(indIm,indIm,:)];
        
        indRe = [1:2*tau:TwoDtau-1];
        indIm = [2:2*tau:TwoDtau];
        covS = [Pfin(indRe,indRe,:),Pfin(indRe,indIm,:); 
                Pfin(indIm,indRe,:),Pfin(indIm,indIm,:)];
        varargout(1) = {covS};
        varargout(2) = {Sfull};
        varargout(3) = {covSfull};
    else
        indRe = [1:2*tau:TwoDtau-1];
        indIm = [2:2*tau:TwoDtau];
        covS = [Pfin(indRe,indRe,:),Pfin(indRe,indIm,:); 
                Pfin(indIm,indRe,:),Pfin(indIm,indIm,:)];
        varargout(1) = {covS};
        varargout(2) = {Sfull};
    end
    
else

    indRe = [1:2*tau:TwoDtau-1];
    indIm = [2:2*tau:TwoDtau];

    S = squeeze(Xfin(1,indRe,:)+i*Xfin(1,indIm,:));
    
    if nargout > 1
        covS = [Pfin(indRe,indRe,:),Pfin(indRe,indIm,:); 
                Pfin(indIm,indRe,:),Pfin(indIm,indIm,:)];
        varargout(1) = {covS};
    end
    
end

fprintf('                                        \r')
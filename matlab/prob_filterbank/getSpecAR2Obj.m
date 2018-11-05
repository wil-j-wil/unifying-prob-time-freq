function [Obj,varargout] = getSpecAR2Obj(theta,vary,bet,specTar,minVar,limCF,limDF);

% function Obj = getSpecAR2Obj(theta,vary,specTar,minVar,limCF,limDF);
% Optional:
% function [Obj,dObj] = getSpecAR2Obj(theta,vary,specTar,minVar,limCF,limDF);
%
% For fitting an AR(2) filter bank to a signal.
% Fits the AR(2) parameters such that the spectra are matched to
% the target spectrum specified in specTar.
%
% The cost function is:
% \sum_t [ log(specAR2_t) + specTar_t/specAR2_t ]
% Which can be derived from the fact that the model is a Gaussian
% with a spectrum parameterised using the AR(2) parameters.
%
% Each AR(2) process is:
% x_{d,t} = \lambda_1 x_{d,t-1} + \lambda_2 x_{d,t-2} + \sigma_d \epsilon
% \epsilon ~ Norm(0,1)
%
% But the function uses an alternative parameterisation in terms of
% the cosine centre frequency, cosine bandwidth and marginal
% variance of the processes as the objective is more simply behaved
% in this space and the constraints on the parameters are simpler
% to enforce. 
%
% INPUTS
% theta = parameters of the AR2FB to be optimised: first D components
%         are the log marginal variances, the next D components are
%         the transformed AR(2) cosine centre frequencies, the final D
%         parameters are the transformed cosine AR(2) bandwidths. Size
%         is therefore [3*D,1]
% vary = white noise variance (if this is set to zero it can result
%        in numerical problems arising from the division
%        spectAR2./specTar) a sensible setting is: vary=max(specTar)*1e-4
% specTar = target spectrum of the data to be fit, size [N,1]
% minVar = minimum of the marginal variances size [D,1]
% limCF = min/max centre-frequencies e.g. [0,1/2], size [D,2]
% limDF = min/max bandwidths e.g. [0,1/2], size [D,2]
% 
% OUTPUTS
% Obj = objective
% OPTIONAL OUTPUT:
% dObj = derivative of the objective wrt the AR(2) parameters, size [3*D,1]
% 
% See Chapter 5 of my thesis (Statistical Models for Natural Sounds
% by R.E.Turner) for more details about AR(2) Filter banks.

D = length(theta)/3;
dVar= exp(theta(1:D));

mVar = minVar+dVar;
cosCF = limCF(:,1)+(limCF(:,2)-limCF(:,1))./(1+exp(-theta(D+1:2*D)));
cosDF = limDF(:,1)+(limDF(:,2)-limDF(:,1))./(1+exp(-theta(2*D+1:3*D)));

N = length(specTar);
Freq = linspace(0,1/2,ceil(N/2));
Freq = [Freq,-Freq([floor(N/2):-1:1])];

% Get the component spectra
a1 = cosDF.^2+4*cosCF.^2;
a2 = -2*cosDF.^2+8*cosCF.^2-2;
z2 = -a1/2-real(sqrt(a1.^2/4-a2+2));

% NB: c1 = (1/2*DF.^2+2*CF.^2+1)*ones(1,N);
oneN = ones(1,N);

c1 = (a1/2+1)*oneN;
c2 = -4*cosCF*oneN;
c3 = (-z2-real(sqrt(z2.^2-4)))*oneN;

lam2 = c3(:,1)*oneN/(-2);

cosCFsq = cosCF.^2*oneN;

preFac = ((lam2-1).^3.*(lam2.^2-1)-16*lam2.^2.*(1+lam2).*cosCFsq) ...
	 ./(2*lam2.*(lam2-1).^3);

cosOm = ones(D,1)*cos(2*pi*Freq);
cos2Om = ones(D,1)*cos(4*pi*Freq);

specsNorm = (c1+c2.*cosOm+cos2Om);

specs = ((mVar*oneN).*preFac)./specsNorm;
spec = sum(specs,1)+vary; % total spectra


Obj1 = sum(log(spec))+sum(specTar'./spec);

Obj2 = bet*sum(mVar);

Obj = (Obj1+Obj2)/N;

%plot(log(specTar),'-k'); hold on; plot(log(spec),'-r')
%keyboard

if nargout>1
  % Derivative of the components wrt parameters

  dCFdtheta = (limCF(:,2)-limCF(:,1))*(1/4)./cosh(theta(D+1:2*D)/2).^2;
  dDFdtheta = (limDF(:,2)-limDF(:,1))*(1/4)./cosh(theta(2*D+1:3*D)/2).^2;
  
  dspecdlogVar = preFac.*(dVar*oneN)./(specsNorm);

  %a1 = cosDF.^2+4*cosCF.^2;
  %a2 = -2*cosDF.^2+8*cosCF.^2-2;

  da1dCF = 8*cosCF;
  da1dDF = 2*cosDF;
  da2dCF = 16*cosCF;
  da2dDF = -4*cosDF;

%  z2 = -a1/2-sqrt(a1.^2/4-a2+2);
  dz2dCF = -da1dCF/2-1/2*(a1.^2/4-a2+2).^(-1/2).*(a1/2.*da1dCF-da2dCF);
  dz2dDF = -da1dDF/2-1/2*(a1.^2/4-a2+2).^(-1/2).*(a1/2.*da1dDF-da2dDF);

%  c3 = (-z2-sqrt(z2.^2-4))*ones(1,N);
  dc3dCF = -dz2dCF-(z2.^2-4).^(-1/2).*z2.*dz2dCF;
  dc3dDF = -dz2dDF-(z2.^2-4).^(-1/2).*z2.*dz2dDF;

% c1 = (a1/2+1)*ones(1,N);
% c2 = -4*cosCF*ones(1,N);
% specsNorm = (c1+c2.*cosOm+cos2Om);
  dspecNormdCF = da1dCF*oneN/2-4*cosOm;
  dspecNormdDF = da1dDF*oneN/2;

%  preFac = ((lam2-1)^3*(lam2^2-1)-16*lam2^2*cosCF^2*(1+lam2)) ...
%	 /(2*lam2*(lam2-1)^3);

dpreFacdlam2 = (3*(lam2-1).^2.*(lam2.^2-1)+2*lam2.*(lam2-1).^3 ...
	    -32*lam2.*(1+lam2).*cosCFsq-16*lam2.^2.*cosCFsq) ...
	    ./(2*lam2.*(lam2-1).^3) ...
	    -preFac.*(2*(lam2-1).^3+6*lam2.*(lam2-1).^2) ...
	    ./(2*lam2.*(lam2-1).^3);
dlam2dc3 = -1/2;

dpreFacdCF = dlam2dc3*dpreFacdlam2.*(dc3dCF*oneN) ...
             -32*lam2.^2.*(cosCF*oneN).*(1+lam2)./(2*lam2.*(lam2-1).^3);

dpreFacdDF = dlam2dc3*dpreFacdlam2.*(dc3dDF*oneN);

% specs = ((Var.*preFac)*ones(1,N))./specsNorm;

dspecdCF = (mVar*ones(1,N)).*dpreFacdCF./specsNorm ...
      - specs./specsNorm.*dspecNormdCF;

dspecdDF = (mVar*ones(1,N)).*dpreFacdDF./specsNorm ...
      - specs./specsNorm.*dspecNormdDF;

dspecdtheta = zeros(3*D,N);
dspecdtheta(1:D,:) = dspecdlogVar;
dspecdtheta(D+(1:D),:) = dspecdCF.*(dCFdtheta*oneN);
dspecdtheta(2*D+(1:D),:) = dspecdDF.*(dDFdtheta*oneN);
  
  dObj1 = sum(ones([3*D,1])*(1./spec.*(1-specTar'./spec)).*dspecdtheta,2);

  dObj2 = [bet*dVar;zeros(2*D,1)];
  
  dObj = [dObj1+dObj2]/N;
  varargout{1}= dObj;

%  Obj = sum(c3(:,1));
%  dObj = [zeros(D,1);dc3dCF;dc3dDF];

%  Obj = sum(z2);
%  dObj = [zeros(D,1);dz2dCF;dz2dDF];
  
%  Obj = sum(specsNorm(:));
%  dObj = [zeros(D,1);sum(dspecNormdCF,2);sum(dspecNormdDF,2)];

% Obj = sum(preFac(:));
% dObj = [zeros(D,1); ...
% 	sum(dpreFacdCF.*(dCFdtheta*oneN),2); ...
% 	sum(dpreFacdDF.*(dDFdtheta*oneN),2)];

end


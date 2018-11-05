function specs = getSpecAR2cosFreq(cosCF,cosDF,mVar,Freq);

% specs = getSpecAR2cosFreq(cosCF,cosDF,mVar,Freq);
%
% Converts from the cosine frequency/marginal variance
% representation into frequency space
%
% INPUTS
% cosCF = cos centre frequencies size [D,1]
% cosDF = cos bandwidths size [D,1]
% mVar = marginal variance of the processes [D,1]
% Freq = frequencies over which to evaluate the spectra [1,N]
% 
% OUTPUTS
% specs = spectra of the processes [D,N]
% 
% See Chapter 5 of my thesis (Statistical Models for Natural Sounds
% by R.E.Turner) for more details about AR(2) Filter banks.

D = length(mVar);
N = length(Freq);

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



function [fmax,df,varMa] = probSpec2freq(om,lamx,varx)

% function [fmax,df,varMa] = probSpec2freq(om,lamx,varx)
%
% Finds the spectral parameters of a probabilistic spectrogram.
% The bandwidth returned (df) is roughly the full width at half
% maximum, however since lamx<~0.17 don't ever fall below half the
% maximum, this is capped to 1/2 for lambdas less than this.
%
% INPUTS
% lamx = dynamical AR parameters [D,1]
% varx = dynamical noise parameters [D,1]
% om = mean frequencies of the sinusoids [D,1]
%
% OUTPUTS
% fmax = centre frequencies, size [D,1]
% df = bandwidths, size [D,1]
% varMa = marginal variances, size [D,1]
%
% If lamx is so low that

% centre frequency
fmax = om/(2*pi);

% bandwidth
df = acos(2-(lamx.^2+1)./(2*lamx))/(2*pi);

% cap on bandwidth
[temp,lamxMax,temp] = freq2probSpec(0,1/2,1);
ind = lamx<lamxMax;
df(ind)=1/2;

% marginal variance
varMa = varx./(1-lamx.^2);


%[Lam,Var] = probSpec2AR2(om,lamx,varx);
%[fmax,df, varMa] = AR22freq(Lam,Var);
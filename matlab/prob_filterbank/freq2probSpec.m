function [om,lamx,varx] = freq2probSpec(fmax,df,varMa)

% function [om,lamx,varx] = freq2probSpec(fmax,df,varMa)
%
% Finds the parameters of a probabilistic spectrogram which has a
% specrum with centre frequency fmax and bandwidth df, as well as
% marginal variance varMa.
%
% INPUTS
% fmax = centre frequencies, size [D,1]
% df = bandwidths, size [D,1]
% varMa = marginal variances, size [D,1]
%
% OUTPUTS
% lamx = dynamical AR parameters [D,1]
% varx = dynamical noise parameters [D,1]
% om = mean frequencies of the sinusoids [D,1]

om = 2*pi*fmax;
%lamx = cos(2*pi*df)-sqrt(cos(2*pi*df).^2+1);

lamx = 2-cos(2*pi*df)-sqrt(cos(2*pi*df).^2-4*cos(2*pi*df)+3);
varx = varMa.*(1-lamx.^2);

%[Lam,Var] = freq2AR2(fmax,df,varMa);
%[om,lamx,varx] = AR22probSpec(Lam,Var);
function [freqs,spec] = getSpecAR2CFDFmVar(cosCF,cosDF,mVar, ...
						    NumFreqs,RngFreqs);

% function [freqs,spec] = getSpecAR2CFDFmVar(cosCF,cosDF,mVar, ...
%
% Computes the spectrum of a single AR(2) process with cosine
% centre-frequency CF, cosine-frequency bandwidth cosDF and
% marginal variance mVar.
%
%
% INPUTS
% cosCF = cosine of the centre frequencies, scalar
% cosdF = cosine FWHM bandwidths, scalar
% mVar = marginal variance, scalar
% NumFreqs = number of frequencies over which to return
%            the spectra, integer scalar
% RngFreqs = range over which to return the spectra, size [2,1]
% 
% OUTPUTS
% freqs = frequencies at which the spectrum is
%         evaluated [1,NumFreqs]  
% spec = spectrum [1,NumFreqs]  
% 
% See Chapter 5 of my thesis (Statistical Models for Natural Sounds
% by R.E.Turner) for more details about AR(2) Filter banks.

D=1;
freqs = [RngFreqs(1):diff(RngFreqs)/(NumFreqs-1):RngFreqs(2)];

% Get the component spectra
a1 = cosDF.^2+4*cosCF.^2;
a2 = -2*cosDF.^2+8*cosCF.^2-2;
z2 = -a1/2-sqrt(a1.^2/4-a2+2);

c1 = (a1/2+1);
c2 = -4*cosCF;
c3 = (-z2-sqrt(z2.^2-4));
lam2 = c3/(-2);

preFac = ((lam2-1)^3*(lam2^2-1)-16*lam2^2*cosCF^2*(1+lam2))...
	 /(2*lam2*(lam2-1)^3);

cosOm = ones(D,1)*cos(2*pi*freqs);
cos2Om = ones(D,1)*cos(4*pi*freqs);

specsNorm = (c1+c2.*cosOm+cos2Om);

spec = preFac*mVar./specsNorm;


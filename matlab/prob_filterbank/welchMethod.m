function [pg,varpg] = welchMethod(y,numFreq,ovLp)

% function [pg,varpg] = welchMethod(y,numFreq,Overlap)
%
% Computes the Welch Periodogram of a signal 'y' by 1) splitting it up
% into sections of length 'numFreq' and with overlap 'ovLp' 2)
% computing the fft of each chunk and 3) then averaging the result.
% Only returns the positive half of the spectrum. So, for example, the
% variance of the signal will be approximately twice the sum of pg
% i.e. mean(y.^2) = 2*sum(pg). The approximation is due to the fact
% that we are computing the variance from windowed versions of the
% signal. 
%
% FOR A POWER SPECTRAL DENSITY ESTIMATE USE:
% psd = 2*pg/numFreq (i.e pg/(delta_f) where delta_f = (1/2)/numFreq) 
%
% In the future I could implement some windows into this method -
% currently I circularlly symmetrise each chunk of the signal via y ->
% [y;y(end-1:-1:2)] and fft this, but I alternatively I could use a
% windowed version of the signal in this process. In other words,
% currently the spectrum is smoothed with a sinc, but another
% window might work better.
%
% INPUTS
% y = signal, size [T,1]
% numFreq = number of frequencies at which to compute the pg
% ovLp = overlap of each chunk (ovLp<numFreq)
%
% OUTPUTS
% pg = (mean) periodogram, size [numFreq,1]
% varpg = variance in the periodogram estimates, size [numFreq,1]
%
% See the wikipedia entry for details of the method:
% http://en.wikipedia.org/wiki/Welch_method

if ovLp>numFreq
  disp('Error in welchMethod.m line 34: overlap > numFreq');
  pg = NaN;
  varpg = NaN;
  return;
end

T = length(y);
Tc = numFreq;

K = floor((T-ovLp)/(Tc-ovLp));

pg = zeros(numFreq,1);
Epg2 = zeros(numFreq,1);

for k=1:K

  startC = 1+(Tc-ovLp)*(k-1);
  endC = startC+Tc-1;
  yCur = y(startC:endC); 

  specCur = abs(fft([yCur;yCur(end-1:-1:2)])).^2;

  specCur = specCur(1:numFreq)/(2*(numFreq-1))^2;
  
  pg = pg+specCur/K;
  Epg2 = Epg2+specCur.^2/K;
  
end

varpg = (Epg2-pg.^2)/K;
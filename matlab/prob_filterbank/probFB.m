function [Z,varargout] = probFB(y,lamx,varx,om,vary,varargin)
  
  % function Z = probFB(y,lamx,varx,om,vary,varargin)
  %
  % Probabilistic Filter Bank (see Turner 2010, Chapter 5 for details)
  %
  % In the standard mode above the FFT is used to compute the AR
  % filter bank. This is fast. When the optional inputs/output below
  % are added/requested, the kalman mode is used which is slower.
  % 
  % NOTE: THAT THE FFT METHOD GIVES A SLIGHTLY DIFFERENT SOLUTION
  % FROM THE KALMAN BASED METHODS (SEE BELOW). SPECIFICALLY, AT THE
  % START AND END OF THE SIGNAL THERE ARE DISCREPANCIES DUE TO THE
  % CIRCULAR BOUNDARY CONDITIONS ASSUMED. IN GENERAL THOUGH, THE
  % TWO METHODS WILL BE EXTREMELY SIMILAR.
  %
  % With optional inputs and outputs:
  % function [Z,covZ] = probSTFT(y,lamx,varx,om,vary,verbose,KF)
  %
  % INPUTS
  % lamx = dynamical AR parameters [D,1]
  % varx = dynamical noise parameters [D,1]
  % om = mean frequencies of the sinusoids [D,1]
  % vary = oberservation noise
  % y = Data, size [T,1] 
  %
  % OPTIONAL INPUTS:
  % verbose = 1 => verbose output
  % KF = 1 => Kalman Filter (rather than smoothing)
  %  
  % OUTPUTS
  % Z = probabilistic filter bank process mean values (these are
  %     complex), size [D,T]
  % OPTIONAL OUTPUTS:
  % covZ = Covariances of these values, [2D,2D,T]
  %
  % I could modify this to return the sufficient statistics and
  % likelihood if requested by the user - see kalman.m
  
T = length(y);

if nargin<=5
  verbose = 0 ;
else
  verbose = varargin{1};
end

if nargin<=6
  KF = 0 ;
else
  KF = varargin{2};
end

if KF==0 & nargout==1 & length(vary)==1

  Z = probFB_FFT(y,lamx,varx,om,vary);
   
else
  
  % Kalman Smoothing
  [lik,Xfin,Pfin] = kalmanSlowFB(lamx,varx,om,vary,y,verbose,KF);

  % Output
  [Z,covZ] = getFBLDSOutput(Xfin,Pfin);
  
end

if nargout==2
  varargout(1) = {covZ};
else
  varargout(1) = {[]};
end
function [Z,varargout] = kernel_ss_probFB(y,A,Q,C,P0,K,vary,tau,varargin)
  
  % Probabilistic Filter Bank using the kernel-based approach in state-space form
  %
  % Will Wilkinson 20181017, based on code by Richard Turner.
  %
  %
  % With optional inputs:
  % function [Z,covZ] = probSTFT(y,lamx,varx,om,vary,verbose,KF)
  %
  % INPUTS
  % A = state transition matrix [2Dtau,2Dtau] (where tau is order of AR model)
  % Q = process noise covariance [2Dtau,2Dtau]
  % C = measurement matrix [2Dtau,1]
  % P0 = initial noise covariance
  % K = dimensionality (number of frequency channels)
  % vary = oberservation noise
  % y = Data, size [T,1] 
  %
  % OPTIONAL INPUTS:
  % verbose = 1 => verbose output
  % KF = 1 => Kalman Filter (rather than smoothing)
  % slow = 1 => Kalman smoothing, 0 => faster implementation, coming soon..
  %  
  % OUTPUTS
  % Z = probabilistic filter bank process mean values (these are
  %     complex), size [D,T]
  %
  % OPTIONAL OUTPUTS:
  % covZ = Covariances of these values, [2D,2D,T]
  % Zfull = probabilistic filter bank process mean values, including higher order terms (these are
  %     complex), size [tauD,T]
  % covZfull = Covariances of these values, [2tauD,2tauD,T]
 
T = length(y);

if nargin<=8
  verbose = 0 ;
else
  verbose = varargin{1};
end

if nargin<=9
  KF = 0 ;
else
  KF = varargin{2};
end

if nargin<=10
  slow = 0 ;
else
  slow = varargin{3};
end

% Kalman Smoothing
% [lik,Xfin,Pfin] = kernel_ss_kalmanSlowFB(A,Q,C,P0,K,vary,y,verbose,KF);
if slow == 1
    [lik,Xfin,Pfin] = kernel_ss_kalmanSlowFB_rewrite(A,Q,C,P0,K,vary,y,verbose,KF);
else
    error('Faster implementation coming soon. Use Kalman smoothing for now')
end

varargout = cell(1,nargout-1);
% Output in correct format
[Z,varargout{:}] = getFBLDSOutput_tau(Xfin,Pfin,tau);

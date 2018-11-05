function [A,Q,C,R,x0,P0] =  ar2LDSParams(Lam,Var,vary);

% function [A,Q,C,R,x0,P0] =  ar2LDSParams(Lam,Var,vary);
%
% Converts the parameters into the form required by the
% Kalman smoothing algorithm, i.e. from AR(2) with D
% hidden variables to AR(1) with 2D hidden variables.
%
% INPUTS
% Lam = dynamical weights, size [D,2]
% Var = innovation noise variance, size [D,1]
% vary = observation noise
%
% OUPUTS
% A = dynamics matrix, size [2*D,2*D]
% Q = dynamics noise, size [2*D,2*D]
% C = emission matrix, size [1,2*D]
% R = emission noise
% 
% see kalman.m for what these new outputs are used for
% see test_ar2LDSParams.m for tests

  D = length(Var);

  A = zeros(2*D);
  Q = zeros(2*D);
  C = repmat([1,0],[1,D]);

  if length(vary)==1
    R = vary;
  else
    R = reshape(vary,[1,1,length(vary)]);
  end
  
  x0 = zeros(2*D,1);
  
  for d=1:D    
    ind = 2*(d-1)+[1:2];
    
    % Dynamics
    A(ind,ind) = [Lam(d,:);1,0];
    Q(ind(1),ind(1)) = Var(d);
  
    % Prior
    autoCor = getAutoCorARTau(Lam(d,:)',Var(d),2);
    P0(ind,ind) = [autoCor(1),autoCor(2);autoCor(2),autoCor(1)]; 
  end
  
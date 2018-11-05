function [Y,Z] = samplePFB(Lam,Var,Om,vary,T);
  
  % [Y,Z] = samplePFB(Lam,Var,Om,vary,T);
  %
  % Samples T time-steps from a complex AR(2) Filter Bank process.
  % x_{1,t,d} = lam_d x_{1,t-1,d} +\eta_{1,t,d} \varx_d^{1/2}
  % x_{2,t,d} = lam_d x_{2,t-1,d} +\eta_{2,t,d} \varx_d^{1/2}
  % \eta_{i,t,d} ~ \Norm(0,1)
  % z_{t,d} =  exp(i om_d)*(x_{1,t,d}+i x_{2,t,d})
  % y_t = \sum_{d} z_{t,d} + \epsilon_t vary^{1/2}
  %
  % Uses the spectrum of the complex AR(2) process and filters white Gaussian
  % noise. See samplePFBslow.m for the auto-regressive version.
  %
  % Inputs
  % Lam = initial dynamical parameters [D,1]
  % Var = initial dynamic noise variance [D,1]
  % Om = centre frequencies [D,1]
  % vary = initial observation noise 
  % T = number of time-steps to sample
  %
  % Outputs
  % Y = observations [T,1] 
  % X = latent variables [T,D]
  %
  % see Turner, 2010 Chapter 5 for details of the complex AR(2) filter bank 
  % see test_samplePFB.m for the unit tests

  [Y,X] = samplePSTFT(Lam,Var,Om,vary,T);

  [T,D] = size(X);
  Z = zeros(T,D);
  
  for d=1:D
    % Generate the observations from the xs 
    Z(:,d) = X(:,d).*exp(i*Om(d)*[1:T]');
  end
      
  
 
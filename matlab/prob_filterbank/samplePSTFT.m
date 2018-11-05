function [Y,X] = samplePSTFT(Lam,Var,Om,vary,T);
  
  % [Y,X] = samplePSTFT(Lam,Var,Om,vary,T);
  %
  % Samples T time-steps from a complex AR(2) Filter Bank process.
  % x_{1,t,d} = lam_d x_{1,t-1,d} +\eta_{1,t,d} \varx_d^{1/2}
  % x_{2,t,d} = lam_d x_{2,t-1,d} +\eta_{2,t,d} \varx_d^{1/2}
  % \eta_{i,t,d} ~ \Norm(0,1)
  % y_t = \sum_{d} exp(i om_d)*(x_{1,t,d}+i x_{2,t,d})
  %
  % Uses the spectrum of the complex AR(2) process and filters white Gaussian
  % noise. See samplePSTFTslow.m for the auto-regressive version.
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
  % see test_samplePSTFT.m for the unit tests
  
  D = length(Var);  
  X = repmat(0,[T,D]);
    
  % Sample from the AR process by filtering coloured noise
  
  RngFreqs = [0,1/2];
  Y = randn(T,1)*sqrt(vary);    
  for d=1:D

    tau = 5*ceil(2*pi/Om(d)); % if tau=0 the first and last samples
                          % will be correlated - the fft induces
                          % circular correlations
			  
    Tx = 2^ceil(log2(T+tau)); % make a power of two so fft is fast
    
    NumFreqs = Tx/2+1;
    
    % Get spectrum centre frequency
    [Freqs,spec] = getCompSpecPFB(Lam(d),Var(d),0,NumFreqs,RngFreqs);     
    spec = [spec,spec(end-1:-1:2)]';
    
%    x1 = ifft(sqrt(spec).*fft(randn(Tx,1)));
%    x2 = ifft(sqrt(spec).*fft(randn(Tx,1)));
%    X(:,d) = x1(1:T) + i*x2(1:T);     

    z =  ifft(sqrt(spec).*fft(randn(Tx,1)+i*randn(Tx,1)));
    X(:,d) = z(1:T);     
  
    % Generate the observations from the xs 
    Y = Y+real(X(:,d).*exp(i*Om(d)*[1:T]'));
  end
      
  
 
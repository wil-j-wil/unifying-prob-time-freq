function [Freqs,Spec] = getCompSpecPFB(lams,varx,om,NumFreqs,RngFreqs)

  % function [Freqs,Spec] =
  % getCompSpecPFB(lams,varx,om,NumFreqs,RngFreqs)
  %
  % Computes the spectra of the following process:
  % x_{1,t} = lam x_{1,t-1} +\eta_{1,t} \varx^{1/2}
  % x_{2,t} = lam x_{2,t-1} +\eta_{2,t} \varx^{1/2}
  % \eta_{i,t} ~ \Norm(0,1)
  % y_t = exp(i om)*(x_{1,t}+i x_{2,t})
  %
  % INPUTS
  % lams = AR(1) dynamical parameter 
  % varx = AR(1) dynamical noise 
  % om = frequencies
  % NumFreqs = number of frequencies over which to return
  %            the spectra
  % RngFreqs = range over which to return the spectra
  %
  % OUTPUTS
  % Freqs = frequencies at which the spectrum is
  %         evaluated [1,NumFreqs]  
  % Spec = spectrum [1,NumFreqs]  
  
Freqs = [RngFreqs(1):diff(RngFreqs)/(NumFreqs-1): ...
         RngFreqs(2)];

Omegas = 2*pi*Freqs;

Spec = varx./(1+lams.^2-2*lams*cos(Omegas-om));


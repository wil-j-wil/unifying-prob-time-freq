function [CF,dF, mVar] = AR22freq(Lam,Var);
  
  % function [CF,dF, mVar] = AR22freq(Lam,Var);
  % 
  % Converts the parameters of a set of AR(2) processes into
  % associated centre frequencies (CF) and bandwidths
  % (dF) 
  %
  % INPUTS
  % Lam = dynamical parameters, size [D,2] (D = number of AR2 processes)
  % Var = innovations variance, size [D,1]
  %
  % OUTPUTS
  % CF = centre frequencies, size [D,1]
  % dF = FWHM bandwidths, size [D,1]
  % mVar = marginal variances, size [D,1] 

  FSamp = 1;
  D = length(Var);
  
  lam1 = Lam(:,1);
  lam2 = Lam(:,2);
  
  cosOm = lam1.*(lam2-1)./(4*lam2); 
  CF = FSamp/(2*pi) * acos(cosOm);
  c1 = -8-16*cosOm.^2;
  c2 = -4./lam2.*(1+16*lam2.^2.*cosOm.^2./(lam2-1).^2+lam2.^2);

  dcosOm = 1/4*sqrt(c1+c2);
  cosOmPlus = cosOm + dcosOm; 
  cosOmMinus = cosOm - dcosOm;
  dF = FSamp*abs(acos(cosOmPlus)-acos(cosOmMinus))/(2*pi);
  
  mVar = getmVarAR2(Lam,Var);
%  mVar = (1-lam2).*Var./(1-lam1.^2-lam2.^2-lam2-lam1.^2.*lam2+lam2.^3);
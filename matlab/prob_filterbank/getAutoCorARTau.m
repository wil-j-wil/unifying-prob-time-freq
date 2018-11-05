function autoCor = getAutoCorARTau(lam,varx,T) 
  
  % function autoCor = getAutoCorARTau(lam,varx,T)
  %
  % Gets the autocorrelation of an AR(Tau) process by
  % solving a system of linear equations
  % 
  % INPUTS
  % lams = AR(tau) dynamical parameter [tau,1]
  % varx = AR(tau) dynamical noise 
  % T = number of time-steps to compute the
  %     auto-correlation over
  %  
  % OUTPUTS
  % autoCor = Autocorrealtion <x_{t}x_{t-t'}> for t' =0
  %           to T [T,1]

tau = length(lam);
LAM = lam*lam';

L = zeros(tau+1);
L(1,:) = [1+lam'*lam,-2*lam'];
L(2:tau+1,2:tau+1) = -eye(tau);

% The top Row
for t=2:tau
  L(1,t) = L(1,t)+2*sum(diag(LAM,t-1));
end

% The upper left triangle of lams
for t = 2:tau+1
  tstart = t-1;
  L(t,1:tau-t+2) = L(t,1:tau-t+2)+lam(t-1:tau)';
end

% The lower right triangle of lams
for t = 3:tau+1
  L(t,2:t-1) = L(t,2:t-1)+lam(t-2:-1:1)';
end

Var = [varx;zeros(tau,1)];

% Solving the system of linear equations by matrix inversion
autoCor = L\Var;


if T<=tau+1
  autoCor = autoCor(1:T,1);
else
  Textra = T-(tau+1);
  autoCor  = [autoCor;zeros(Textra,1)];%
  for t = 1:Textra
    autoCor(t+tau+1,1) = lam'*autoCor(t+tau:-1:t+1,1);
  end 
end
  
  
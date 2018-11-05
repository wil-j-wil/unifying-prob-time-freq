function [A,Q,H,Pinf,K,tau1] = get_disc_model(lamx,varx,omega,D,kernel,se_approx_order)



  % Define the hyperparameters
  dt = 1;  % step size is 1 sample, regardless of sample rate
  w = omega;  % om
  
  if strcmp(kernel,'exp')
    lengthScale = 1 ./ lamx;
  elseif strcmp(kernel,'matern32')
    lengthScale = sqrt(3) ./ lamx;
  elseif strcmp(kernel,'matern52')
    lengthScale = sqrt(5) ./ lamx;
  end
  magnSigma2 = varx;
  
  
  
  Qc1 = zeros(D,1);
  F1=[];L1=[];H1=[];Pinf1=[];
  for d=1:D
    % Construct the continuous-time SDE (from Solin+Sarkka AISTATS 2014)
    cf_to_ss = str2func(strcat('cf_',kernel,'_to_ss'));
    [F1d,L1d,Qc1d,H1d,Pinf1d] = cf_to_ss(magnSigma2(d), lengthScale(d), se_approx_order);
    F1 = blkdiag(F1, F1d);
    L1 = vertcat(L1,L1d);
    Qc1(d) = Qc1d;
    H1 = horzcat(H1,H1d);
    Pinf1 = blkdiag(Pinf1,Pinf1d);
  end
  tau1 = length(L1d); % tau = model order (1 for Exponential, 2 for Matern 3/2, etc.)
  
    
  % Construct the continuous-time SDE: periodic (just a sinusoid) (from Solin+Sarkka AISTATS 2014)
  tau2=2; % real + imaginary
  F2=[];L2=[];Qc2=[];H2=[];
  for d=1:D
    F2 = blkdiag(F2,[0 -w(d); w(d) 0]);
    L2 = blkdiag(L2,eye(tau2));
    Qc2 = blkdiag(Qc2,zeros(tau2));
    H2 = horzcat(H2,[1 0]);
  end


%%
  % The product of the two (from Solin+Sarkka AISTATS 2014)
%   F = kron(F1,eye(tau2)) + kron(eye(tau1),F2);  % <- first order approach
%   L = kron(L1,eye(tau2));
%   Qc = kron(Qc1,eye(tau2));
  F2_kron=[];L=[];Qc=[];
  for d=1:D  % for higher-order models we must iterate to stack the kronecker products along the diagonal
    idx1 = tau1*(d-1)+1:tau1*d;
    idx2 = tau2*(d-1)+1:tau2*d;
    F2d = F2(idx2,idx2);
    F2d_kron = kron(eye(tau1),F2d);
    F2_kron = blkdiag(F2_kron,F2d_kron);
    L = blkdiag(L, kron(L1(idx1),L2(idx2,idx2)));
    Qc = blkdiag(Qc, kron(Qc1(d),L2(idx2,idx2)));
  end
  F = kron(F1,eye(tau2)) + F2_kron;
  H = kron(H1,[1 0]);
  Pinf = kron(Pinf1,eye(tau2));
  
  K = D*tau1; % model order
  
  % Solve the discrete-time state-space model (the function is from the EKF/UKF toolbox)
  [A,Q] = lti_disc(F,L,Qc,dt);
  
end
% example script walking through the steps for converting a spectral
% mixture GP model to state space form.

clear;

% filterbank code
addpath('prob_filterbank');

% Specify where to load the data from
soundPath = '../audio/speech/';

filename = 'speech0_female';

% number of frequency channels
D = 20;

% sample rate
fs = 16000;


% Choose a kernel
kernel = 'matern32';  % 'exp', 'matern32', 'matern52', 'se'
se_approx_order = 4;


% measurement noise
vary = 1e-8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load signal and pre-process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[y,fs_] = audioread([soundPath,filename,'.wav']); % reads in the file
RngLimTest = round([1, fs_]); % 1 second of audio
yTest = y(RngLimTest(1):RngLimTest(2),1); 
yTest = resample(yTest, fs, fs_); % downsample the input
sig_var = sqrt(var(yTest));
yTest = yTest/sig_var; % rescale the input to unit variance
T = length(yTest);

% for now don't fit the filterbank to the spectrum, just use log-spacing
dfFrac = 1/20; % best results dfFrac = 1/20
  fmax = logspace(log10(1/200),log10(0.3),D)';
  % [omega,psi,ro] = [instantaneous frequency, shrinkage/model variance, noise variance]
  [omega,lamx,varx] = freq2probSpec(fmax,fmax*dfFrac,ones(D,1)/D);  % lamx=psi, varx=ro
  
% standard model spectrogram computed here
[ZStandard, ~] = probFB(yTest,lamx,varx,omega,vary); 
AStandard = abs(ZStandard').^2;
yStandard = sum(real(ZStandard'),2);

figure(1); clf
  subplot(2,1,1)
  imagesc(log(AStandard)')
  set(gca,'YDir','normal')
  title('spectrogram of y (standard model)')
  subplot(2,1,2)
  plot(yStandard)
  title('signal y (standard model)')

%%

  % Define the hyperparameters
  dt = 1;  % step size is 1 sample, regardless of sample rate
  w = omega;  % om
  
  if strcmp(kernel,'exp')
    lengthScale = -dt ./ log(lamx);
  elseif strcmp(kernel,'matern32')
    lengthScale = -dt * sqrt(3) ./ log(lamx);
  elseif strcmp(kernel,'matern52')
    lengthScale = -dt * sqrt(5) ./ log(lamx);
  end
  magnSigma2 = varx ./ (1 - exp(-2 * dt ./ lengthScale));
  
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
  
  
  % Solve the discrete-time state-space model (the function is from the EKF/UKF toolbox)
  % Note that the matrix exponential here gives a slightly different A to the original
  % method (which uses the rotation metrix directly), but results are very similar.
  [A,Q] = lti_disc(F,L,Qc,dt);
  

  
%%
  % kernel-based model spectrogram computed here
  verbose = 1;
  filter_only = 0;
  % Kalman Smoothing
  % Requesting all the outputs is memory intensive. Z only contains the (D) observable
  % dimensions. Zfull contains the higher-order terms too.
%   ZKernel = kernel_ss_probFB(yTest,A,Q,H,Pinf,D*tau1,vary,tau1,verbose,filter_only);
  [ZKernel, covZ, ZKernelfull, covZfull] = kernel_ss_probFB(yTest,A,Q,H,Pinf,D*tau1,vary,tau1,verbose,filter_only,1);
  AKernel = abs(ZKernel').^2;
  yKernel = sum(real(ZKernel'),2);

  figure(2); clf
    subplot(2,1,1)
    imagesc(log(AKernel)')
    set(gca,'YDir','normal')
    title(sprintf('spectrogram of y (%s model)', kernel))
    subplot(2,1,2)
    plot(yKernel)
    title(sprintf('signal y (%s model)', kernel))
    
    

%%
 % Compare the energy in each frequency channel between standard model and
 % the new one. Notice how higher-order models => narrow bands.
    freqs = omega*fs/(2*pi);
    figure(3); clf
    for d=1:D
    fprintf('channel %d, centre freq.: %6.1f Hz \n', d, freqs(d))
    channely = sum(real(ZStandard(d,:)'),2);
    ffty = fft(channely)./numel(channely);
    channelykern = sum(real(ZKernel(d,:)'),2);
    fftykern = fft(channelykern)./numel(channelykern);
    subplot(2,1,1) 
    plot(real(ffty(1:T/2)))
    ylim([-0.04, 0.04])
    xlim([0, freqs(end)])
    hold on
    title(sprintf('Standard first order model - energy in freq. channels 1-%d', d))
    subplot(2,1,2)
    plot(real(fftykern(1:T/2)))
    ylim([-0.04, 0.04])
    xlim([0, freqs(end)])
    hold on
    title(sprintf('%s model - energy in freq. channels 1-%d', kernel, d))
    xlabel('Frequency (Hz)')
    pause
    end
    fprintf('\n')
 
 
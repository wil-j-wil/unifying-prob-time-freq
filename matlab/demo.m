% example script for fitting a prob. time frequency model in the frequency
% domain. Here we pre-optimise using the exponential kernel, which gives
% much better results.

clear; close all;

% filterbank code
addpath('prob_filterbank/');

% Specify where to load the data from
soundPath = '../audio/speech/';

% load signal
File = 'speech0_female'; % Name of file to load
fs_ = 16000; % sampling rate of file

% DS = 1; % down sample further if requested
D = 40; % number channels (don't set too high)

kernel2 = 'matern32';
se_approx_order = 4;




L = 1; % number of gap lengths to consider

gapLim = [10,320];
% y_max_length = 1; % seconds
% gapPos = [2000,4000,6000,8000,10000,12000,14000]; % manually set gap positions
gapPos = [1500,5000,7000,9000,13000,18000]; % set manually to non-silent regions
% gapPos(gapPos>=y_max_length*fs_-gapLim(end))=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load signal and pre-process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[y,fs] = audioread([soundPath,File,'.wav']); % reads in the file
% yTest = y(1:y_max_length*fs,1); 
yTest = resample(y, fs_, fs); % downsample the input
yTest = yTest(1:min(18300,length(yTest)));  % more than 19000 takes up too much memory if Kalman smoother is used
fs = fs_;
normaliser = sqrt(var(yTest));
yTest = yTest/normaliser; % rescale the input to unit variance
T = length(yTest);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pre-optimise the filter bank with exponential kernel (probabilistic phase vocoder)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yTrain = yTest;

%if trainFilter==1
  % Learn properties of the filters (centre frequency and width)
  opts = struct;
  opts.verbose = 1; % view plots of the fitting process
  opts.minT = 100;
  opts.maxT = 1000;
  opts.numIts = 15;
  opts.numLevels = 30;
  opts.bet = 750;  % increase if filters bunch together
  opts.reassign = 0;  % set to 1 if filters really bunch together
  
  [Var1Fit,Lam1Fit,om1Fit,InfoFit1] = fit_probSTFT_SD(yTrain,D,'exp',opts); % trains filters to match the spectrum
  
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Then tune parameters with new kernel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if trainFilter==1
  % Learn properties of the filters (centre frequency and width)
  opts = struct;
  opts.verbose = 1; % view plots of the fitting process
  opts.minT = 500;
  opts.maxT = 1000;
  opts.numIts = 15;
  opts.numLevels = 20;
  opts.bet = 750;  % increase if filters bunch together
  opts.reassign = 0;  % set to 1 if filters really bunch together
  
  opts.theta_init = [Var1Fit;Lam1Fit;om1Fit];
  [Var2Fit,Lam2Fit,om2Fit,InfoFit2] = fit_probSTFT_SD(yTrain,D,kernel2,opts); % trains filters to match the spectrum

  
 %%
  % plots the sprectrum of the filter
  figH1 = plot_pSTFT_kern_cts(Var2Fit,om2Fit,Lam2Fit,kernel2,fs,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute the spectrogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$$$$

[om2Fit, om_ind] = sort(om2Fit);
Lam2Fit = Lam2Fit(om_ind);
Var2Fit = Var2Fit(om_ind);
[A2,Q2,H2,Pinf2,K2,tau2] = get_disc_model(Lam2Fit,Var2Fit,om2Fit,D,kernel2,se_approx_order);

% clean spectrogram computed here
disp('computing signal spectrogram')
ZTest = kernel_ss_probFB(yTest,A2,Q2,H2,Pinf2,K2,0,tau2,0,0,1);
ATest = abs(ZTest').^2;
yTest = sum(real(ZTest'),2);

  figure(2); clf
    subplot(2,1,1)
    imagesc(log(ATest)')
    set(gca,'YDir','normal')
    title(sprintf('spectrogram of y (%s model)', kernel2))
    subplot(2,1,2)
    plot(yTest)
    title(sprintf('signal y (%s model)', kernel2))
 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

slow = 1; % faster implementation coming soon

NGaps = length(gapPos);

gaps = ceil(linspace(gapLim(1),gapLim(2),L));

numMethods = 2; % includes the un-denoised signal

snr_a = zeros(L,D,numMethods);
snr_loga = zeros(L,D,numMethods);
snr_y = zeros(L,numMethods);
pesq_y = zeros(L,numMethods);

% Ys = zeros(T,numMethods,L);
cnt = 0;

for l=1:L
  
  % display progress to user
  cnt= cnt+1;
  disp(['Progress ',num2str(cnt),'/',num2str(L)]);

  % set the state 
  randn('state',1);
  
  % signal with missing data
  yGap = yTest;
  ind = [];
  for ng=1:NGaps
    ind = [ind,gapPos(ng)+[-ceil(gaps(l)/2):+ceil(gaps(l)/2)]];
  end
  
  yGap(ind) = NaN;
  varys = 0;
  if slow == 1
    yGap(ind) = 0;
    varys = 1e-4*ones(T,1);
    varys(ind) = 1e5;
  end
  
  disp('reconstructing signal')
  ZRecon = kernel_ss_probFB(yGap,A2,Q2,H2,Pinf2,K2,varys,tau2,0,0,slow);
  yRecon = sum(real(ZRecon'),2);
  
  %%% UNCOMMENT THIS SECTION TO PERFORM SPECTROGRAM RECONSTRUCTION %%%
  %{
  disp('computing spectrogram with gaps')
  ZGap = kernel_ss_probFB(yGap,A2,Q2,H2,Pinf2,K2,0,tau2,0,0,slow);
  AGap = abs(ZGap').^2;
  
  disp('reconstructing spectrogram')
  ZTemp = kernel_ss_probFB(yRecon,A2,Q2,H2,Pinf2,K2,0,tau2,0,0,slow);
  ARecon = abs(ZTemp').^2;
  
  % figure out how well we have reconstructed the spectrogram
  
  snr_a(l,:,1) = snr(ATest(ind,:),AGap(ind,:));
  snr_a(l,:,2) = snr(ATest(ind,:),ARecon(ind,:)); 
  
  snr_loga(l,:,1) = snr(log10(loud_floor(ATest(ind,:),deltaSNR)), ...
			log10(loud_floor(AGap(ind,:),deltaSNR)));
  snr_loga(l,:,2) = snr(log10(loud_floor(ATest(ind,:),deltaSNR)), ...
			log10(loud_floor(ARecon(ind,:),deltaSNR)));
  
  snr_loga_mn = squeeze(mean(snr_loga,2))
  snr_a_mn = squeeze(mean(snr_a,2))
  %}
  %%%                                                              %%%
  
%   Ys(:,:,l) = single([yGap,yRecon]);  

  deltaSNR = 1e-5; % apply loudness floor

  snr_y(l,1) = snr(yTest(ind,:),yGap(ind,:));
  snr_y(l,2) = snr(yTest(ind,:),yRecon(ind,:));

%   pesq_y(l,1) = pesq(yTest, yGap, fs_);
%   pesq_y(l,2) = pesq(yTest1, yRecon1, fs_);

  snr_y
%   pesq_y
  

end


%
grey = [0.8 0.8 0.8];
figure(3); clf
subplot(2,1,1)
plot(yTest,'Color','r')
hold on
plot(yGap,'b')
title('actual signal')
y_limits = ylim;
subplot(2,1,2)
plot(yRecon,'b')
title(sprintf('reconstruction with %s model', kernel2))
% title(sprintf('reconstruction with %s model', 'exp'))
ylim(y_limits)

t = linspace(1,T,T);
ind1gap = gapPos(1)+[-ceil(gaps(L)/2):+ceil(gaps(L)/2)];
ind2gap = gapPos(2)+[-ceil(gaps(L)/2):+ceil(gaps(L)/2)];
ind1 = gapPos(1)+[-2.5*ceil(gaps(L)/2):+2.5*ceil(gaps(L)/2)];
ind2 = gapPos(2)+[-2.5*ceil(gaps(L)/2):+2.5*ceil(gaps(L)/2)];
figure(4); clf
subplot(2,2,1)
plot(t(ind1),yTest(ind1),'Color',grey)
hold on
plot(t(ind1gap),yTest(ind1gap), 'k')
title('actual signal')
y_limits = ylim;
subplot(2,2,3)
plot(t(ind1),yRecon(ind1),'Color',grey)
hold on
plot(t(ind1gap),yRecon(ind1gap), 'b')
title(sprintf('reconstruction with %s model', kernel2))
% title(sprintf('reconstruction with %s model', 'exp'))
ylim(y_limits)
subplot(2,2,2)
plot(t(ind2),yTest(ind2),'Color',grey)
hold on
plot(t(ind2gap),yTest(ind2gap), 'k')
title('actual signal')
y_limits = ylim;
subplot(2,2,4)
plot(t(ind2),yRecon(ind2),'Color',grey)
hold on
plot(t(ind2gap),yRecon(ind2gap), 'b')
title(sprintf('reconstruction with %s model', kernel2))
% title(sprintf('reconstruction with %s model', 'exp'))
ylim(y_limits)

%%
% sound(yRecon*normaliser,fs)
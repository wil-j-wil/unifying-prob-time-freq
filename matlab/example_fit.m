% example script for fitting a prob. time frequency model in the frequency
% domain. If using a kernel other than the exponential, then see
% example_fit_pre_opt.m which has a pre-optimisation stage and gives better
% results.

clear; close all;

% filterbank code
addpath('prob_filterbank/');

% Specify where to load the data from
soundPath = '../audio/speech/';

% load signal
File = 'speech0_female'; % Name of file to load
fs_ = 22050; % sampling rate of file

% DS = 1; % down sample further if requested
D = 30; % number channels (don't set too high)

kernel = 'exp';
se_approx_order = 4;





L = 1; % number of gap lengths to consider

gapLim = [10,500];
y_max_length = 1; % seconds
gapPos = [2000,4000,6000,8000,10000,12000,14000]; % manually set gap positions
gapPos(gapPos>=y_max_length*fs_-gapLim(end))=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load signal and pre-process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[y,fs] = audioread([soundPath,File,'.wav']); % reads in the file
yTest = y(1:y_max_length*fs,1); 
yTest = resample(yTest, fs_, fs); % downsample the input
fs = fs_;
normaliser = sqrt(var(yTest));
yTest = yTest/normaliser; % rescale the input to unit variance
T = length(yTest);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Train the filter bank
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yTrain = yTest;

%if trainFilter==1
  % Learn properties of the filters (centre frequency and width)
  opts.verbose = 1; % view plots of the fitting process
  opts.minT = 100;
  opts.maxT = 1000;
  opts.numIts = 10;
  opts.numLevels = 30;
  opts.bet = 750;  % increase if filters bunch together
  opts.reassign = 0;  % set to 1 if filters really bunch together
  
  [Var1Fit,Lam1Fit,omFit,InfoFit] = fit_probSTFT_SD(yTrain,D,kernel,opts); % trains filters to match the spectrum

  
 %%
  % plots the sprectrum of the filter
%     figH1 = plot_pSTFT(Var1Fit,omFit,Lam1Fit,fs,1);
    figH1 = plot_pSTFT_kern_cts(Var1Fit,omFit,Lam1Fit,kernel,fs,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute the spectrogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$$$$

[A1,Q1,H1,Pinf1,K1,tau1] = get_disc_model(Lam1Fit,Var1Fit,omFit,D,kernel,se_approx_order);

% clean spectrogram computed here
% ZTest = kernel_ss_probFB(yTest,Lam1,Var1,om,0);
disp('computing signal spectrogram')
ZTest1 = kernel_ss_probFB(yTest,A1,Q1,H1,Pinf1,K1,0,tau1,0,0,1);
ATest1 = abs(ZTest1').^2;
yTest1 = sum(real(ZTest1'),2);

  figure(2); clf
    subplot(2,1,1)
    imagesc(log(ATest1)')
    set(gca,'YDir','normal')
    title(sprintf('spectrogram of y (%s model)', kernel))
    subplot(2,1,2)
    plot(yTest1)
    title(sprintf('signal y (%s model)', kernel))
 


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

Ys = zeros(T,numMethods,L);
cnt = 0;

for l=1:L
  
  % display progress to user
  cnt= cnt+1;
  disp(['Progress ',num2str(cnt),'/',num2str(L)]);

  % set the state 
  randn('state',1);
  
  % signal with missing data
  yGap = yTest;
  yGap1 = yTest1;
  ind = [];
  for ng=1:NGaps
    ind = [ind,gapPos(ng)+[-ceil(gaps(l)/2):+ceil(gaps(l)/2)]];
  end
  
  yGap(ind) = NaN;
  yGap1(ind) = NaN;
  varys = 0;
  if slow == 1
    yGap(ind) = 0;
    yGap1(ind) = 0;
    varys = 1e-4*ones(T,1);
    varys(ind) = 1e5;
  end
  
  disp('computing spectrogram with gaps')
%   ZGap = kernel_ss_probFB(yGap,Lam1,Var1,om,0); 
  ZGap1 = kernel_ss_probFB(yGap1,A1,Q1,H1,Pinf1,K1,0,tau1,0,0,slow);
  AGap1 = abs(ZGap1').^2;
  
  
  
  % denoise using model 1
  disp('reconstructing spectrogram with model 1')
  ZRecon1 = kernel_ss_probFB(yGap1,A1,Q1,H1,Pinf1,K1,varys,tau1,0,0,slow);
%   ZRecon1 = kernel_ss_probFB(yGap1,A1,Q1,H1,Pinf1,K1,0,tau1);
  yRecon1 = sum(real(ZRecon1'),2);
%   ZTemp = kernel_ss_probFB(yReconGTF_UT,Lam1,Var1,om,0);
  ZTemp = kernel_ss_probFB(yRecon1,A1,Q1,H1,Pinf1,K1,0,tau1,0,0,slow);
  ARecon1 = abs(ZTemp').^2;
  
  Ys(:,:,l) = single([yGap,yRecon1]);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % figure out the amount we have denoised the spectrogram by
  
  snr_a(l,:,1) = snr(ATest1(ind,:),AGap1(ind,:));
  snr_a(l,:,2) = snr(ATest1(ind,:),ARecon1(ind,:));   

  deltaSNR = 1e-5; % apply loudness floor

  % code to check that the floor is only getting rid of a
  % relatively small proportion of the spectrogram
  %
  % figure
  % subplot(2,1,1)
  % imagesc(log(ATest)')
  
  % subplot(2,1,2)
  % [ATest_floor,frac] = loud_floor(ATest,deltaSNR);
  % frac
  % imagesc(log(ATest_floor)')

  
  snr_loga(l,:,1) = snr(log10(loud_floor(ATest1(ind,:),deltaSNR)), ...
			log10(loud_floor(AGap1(ind,:),deltaSNR)));
  snr_loga(l,:,2) = snr(log10(loud_floor(ATest1(ind,:),deltaSNR)), ...
			log10(loud_floor(ARecon1(ind,:),deltaSNR)));

  snr_y(l,1) = snr(yTest(ind,:),yGap(ind,:));
  snr_y(l,2) = snr(yTest1(ind,:),yRecon1(ind,:));

%   pesq_y(l,1) = pesq(yTest, yGap, fs_);
%   pesq_y(l,2) = pesq(yTest1, yRecon1, fs_);

   
  snr_a_mn = squeeze(mean(snr_a,2))
  snr_loga_mn = squeeze(mean(snr_loga,2))
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
subplot(2,1,2)
plot(yRecon1,'b')
title(sprintf('reconstruction with %s model', kernel))

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
subplot(2,2,3)
plot(t(ind1),yRecon1(ind1),'Color',grey)
hold on
plot(t(ind1gap),yRecon1(ind1gap), 'b')
title(sprintf('reconstruction with %s model', kernel))
subplot(2,2,2)
plot(t(ind2),yTest(ind2),'Color',grey)
hold on
plot(t(ind2gap),yTest(ind2gap), 'k')
title('actual signal')
subplot(2,2,4)
plot(t(ind2),yRecon1(ind2),'Color',grey)
hold on
plot(t(ind2gap),yRecon1(ind2gap), 'b')
title(sprintf('reconstruction with %s model', kernel))


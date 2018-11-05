function [snr_y, pesq_y, varargout] = run_missing_data_experiment(File,fs_,D,kernel,L,gapLim,gapPos,slow,verbose)

se_approx_order = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load signal and pre-process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[y,fs] = audioread(File); % reads in the file
yTest = resample(y, fs_, fs); % downsample the input
yTest = yTest(1:min(18300,length(yTest)));  % more than 19000 takes up too much memory is Kalman smoother is used
fs = fs_;
normaliser = sqrt(var(yTest));
yTest = yTest/normaliser; % rescale the input to unit variance
T = length(yTest);
% gapPos(gapPos>=length(y)-gapLim(end))=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pre-optimise the filter bank with exponential kernel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yTrain = yTest;

%if trainFilter==1
  % Learn properties of the filters (centre frequency and width)
  opts = struct;
  opts.verbose = verbose; % view plots of the fitting process
  opts.minT = 100;
  opts.maxT = 1000;
  opts.numIts = 15;
  opts.numLevels = 30;
  opts.bet = 750;  % increase if filters bunch together
  opts.reassign = 0;  % set to 1 if filters really bunch together
  
  [Var1Fit,Lam1Fit,omFit,~] = fit_probSTFT_SD(yTrain,D,'exp',opts); % trains filters to match the spectrum

  
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Then tune parameters with new kernel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if trainFilter==1
  % Learn properties of the filters (centre frequency and width)
  opts = struct;
  opts.verbose = verbose; % view plots of the fitting process
  opts.minT = 500; % don't smooth as much because initialised to previous run
  opts.maxT = 1000;
  opts.numIts = 15;
  opts.numLevels = 20;
  opts.bet = 750;  % increase if filters bunch together
  opts.reassign = 0;  % set to 1 if filters really bunch together
  
  opts.theta_init = [Var1Fit;Lam1Fit;omFit];
  [Var2Fit,Lam2Fit,om2Fit,~] = fit_probSTFT_SD(yTrain,D,kernel,opts); % trains filters to match the spectrum



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute the spectrogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[A1,Q1,H1,Pinf1,K1,tau1] = get_disc_model(Lam1Fit,Var1Fit,omFit,D,'exp',se_approx_order);
[A2,Q2,H2,Pinf2,K2,tau2] = get_disc_model(Lam2Fit,Var2Fit,om2Fit,D,kernel,se_approx_order);

% clean spectrogram computed here
% ZTest = kernel_ss_probFB(yTest,Lam1,Var1,om,0);
disp('computing signal spectrogram')
ZTest1 = kernel_ss_probFB(yTest,A1,Q1,H1,Pinf1,K1,0,tau1,0,0,slow);
yTest1 = sum(real(ZTest1'),2);
ZTest2 = kernel_ss_probFB(yTest,A2,Q2,H2,Pinf2,K2,0,tau2,0,0,slow);
yTest2 = sum(real(ZTest2'),2);
 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% slow = 0;

% L = 5;
NGaps = length(gapPos);

%gaps = ceil(logspace(log10(gapLim(1)),log10(gapLim(2)),L));
gaps = ceil(linspace(gapLim(1),gapLim(2),L));
%gapPos = gapLim(2)+ceil(rand(NGaps,1)*(T-2*gapLim(2)));

numMethods = 3; % includes the un-denoised signal

snr_y = zeros(L,numMethods);
pesq_y = zeros(L,numMethods);


Ys = zeros(T,numMethods,L);
cnt = 0;

% for l=[10,1,5,2,6,9,3,8,4,7]
for l=1:L
  
  % display progress to user
  cnt= cnt+1;
  disp(['Progress ',num2str(cnt),'/',num2str(L)]);

  % set the state 
  randn('state',1);
  
  % signal with missing data
  yGap = yTest;
  yGap1 = yTest1;
  yGap2 = yTest2;
  ind = [];
  for ng=1:NGaps
    ind = [ind,gapPos(ng)+[-ceil(gaps(l)/2):+ceil(gaps(l)/2)]];
  end
  
  yGap(ind) = NaN;
  yGap1(ind) = NaN;
  yGap2(ind) = NaN;
  varys = 0;
  if slow == 1
    yGap(ind) = 0;
    yGap1(ind) = 0;
    yGap2(ind) = 0;
    varys = 1e-4*ones(T,1);
    varys(ind) = 1e5;
  end
  
  % denoise using model 1
  disp('reconstructing spectrogram with model 1')
  ZRecon1 = kernel_ss_probFB(yGap1,A1,Q1,H1,Pinf1,K1,varys,tau1,0,0,slow);
  yRecon1 = sum(real(ZRecon1'),2);
  
  % denoise using model 2
  disp('reconstructing spectrogram with model 2')
  ZRecon2 = kernel_ss_probFB(yGap2,A2,Q2,H2,Pinf2,K2,varys,tau2,0,0,slow);
  yRecon2 = sum(real(ZRecon2'),2);
  
  Ys(:,:,l) = single([yGap,yRecon1,yRecon2]);
  
  yGap(isnan(yGap)) = 1e-5;
%   snr_y(l,1) = snr(yTest(ind,:),yGap(ind,:));
%   snr_y(l,2) = snr(yTest1(ind,:),yRecon1(ind,:));
%   snr_y(l,3) = snr(yTest2(ind,:),yRecon2(ind,:));
  snr_y(l,1) = snr(yTest(ind,:),yTest(ind,:)-yGap(ind,:));
  snr_y(l,2) = snr(yTest1(ind,:),yTest1(ind,:)-yRecon1(ind,:));
  snr_y(l,3) = snr(yTest2(ind,:),yTest2(ind,:)-yRecon2(ind,:));
  
  try
    pesq_y(l,1) = pesq(yTest, yGap, fs_/2);
    pesq_y(l,2) = pesq(yTest1, yRecon1, fs_/2);
    pesq_y(l,3) = pesq(yTest2, yRecon2, fs_/2);
  catch
      disp('PESQ failed')
  end

  snr_y
  pesq_y
  

end

%%
if verbose

    figure(2); clf
    subplot(2,1,1)
    hold on
    title('waveform')
    plot(gaps*1000/fs,pesq_y(:,2),'k--')
    plot(gaps*1000/fs,pesq_y(:,3),'k-')
    legend('exp',kernel)
    xlabel('gap /ms')
    ylabel('PSEQ score')

    subplot(2,1,2)
    hold on
    plot(gaps*1000/fs,snr_y(:,2),'k--')
    plot(gaps*1000/fs,snr_y(:,3),'k-')
    legend('exp',kernel)
    xlabel('gap /ms')
    ylabel('SNR')


    %
    grey = [0.8 0.8 0.8];
    figure(3); clf
    subplot(3,1,1)
    plot(yTest,'Color','r')
    hold on
    plot(yGap,'b')
    title('actual signal')
    y_limits = ylim;
    subplot(3,1,2)
    plot(yRecon1,'b')
    title(sprintf('reconstruction with %s model', 'exp'))
    ylim(y_limits)
    subplot(3,1,3)
    plot(yRecon2,'b')
    title(sprintf('reconstruction with %s model', kernel))
    ylim(y_limits)


    t = linspace(1,T,T);
    ind1gap = gapPos(1)+[-ceil(gaps(L)/2):+ceil(gaps(L)/2)];
    ind2gap = gapPos(2)+[-ceil(gaps(L)/2):+ceil(gaps(L)/2)];
    ind1 = gapPos(1)+[-2.5*ceil(gaps(L)/2):+2.5*ceil(gaps(L)/2)];
    ind2 = gapPos(2)+[-2.5*ceil(gaps(L)/2):+2.5*ceil(gaps(L)/2)];
    figure(4); clf
    subplot(3,2,1)
    plot(t(ind1),yTest(ind1),'Color',grey)
    hold on
    plot(t(ind1gap),yTest(ind1gap), 'k')
    title('actual signal')
    y_limits = ylim;
    subplot(3,2,3)
    plot(t(ind1),yRecon1(ind1),'Color',grey)
    hold on
    plot(t(ind1gap),yRecon1(ind1gap), 'b')
    title(sprintf('reconstruction with %s model', 'exp'))
    ylim(y_limits)
    subplot(3,2,5)
    plot(t(ind1),yRecon2(ind1),'Color',grey)
    hold on
    plot(t(ind1gap),yRecon2(ind1gap), 'r')
    title(sprintf('reconstruction with %s model', kernel))
    ylim(y_limits)
    subplot(3,2,2)
    plot(t(ind2),yTest(ind2),'Color',grey)
    hold on
    plot(t(ind2gap),yTest(ind2gap), 'k')
    title('actual signal')
    y_limits = ylim;
    subplot(3,2,4)
    plot(t(ind2),yRecon1(ind2),'Color',grey)
    hold on
    plot(t(ind2gap),yRecon1(ind2gap), 'b')
    title(sprintf('reconstruction with %s model', 'exp'))
    ylim(y_limits)
    subplot(3,2,6)
    plot(t(ind2),yRecon2(ind2),'Color',grey)
    hold on
    plot(t(ind2gap),yRecon2(ind2gap), 'r')
    title(sprintf('reconstruction with %s model', kernel))
    ylim(y_limits)

    %
    figure(5); clf
    subplot(2,1,1)
    hold on
    title('waveform')
    plot(gaps*1000/fs,pesq_y(:,2)-pesq_y(:,1),'--k')
    plot(gaps*1000/fs,pesq_y(:,3)-pesq_y(:,1),'-k')
    legend('exp',kernel)
    xlabel('gap /ms')
    ylabel('PSEQ improvement')

    subplot(2,1,2)
    hold on
    plot(gaps*1000/fs,snr_y(:,2)-snr_y(:,1),'--k')
    plot(gaps*1000/fs,snr_y(:,3)-snr_y(:,1),'-k')
    legend('exp',kernel)
    xlabel('gap /ms')
    ylabel('SNR improvement')

end

if nargout > 2
    varargout{1} = yTest;
    varargout{2} = yRecon1;
    varargout{3} = yRecon2;
end

end
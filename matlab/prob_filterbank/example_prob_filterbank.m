clear;

% Specify where to load the data from
soundPath = '../../audio/';

% Specify where to save the data to
% saveDir = '/Users/williamwilkinson/Documents/PhD/inProgress/prob_time_freq/data/';

filename = 'stim312_wind';

D = 40;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TESTING
gapPos = 2000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load signal and pre-process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[y,fs] = audioread([soundPath,filename,'.wav']); % reads in the file
RngLimTest = round([1, length(y)]);
yTest = y(RngLimTest(1):RngLimTest(2),1); 
yTest = resample(yTest, 16000, fs); % downsample the input
fs = 16000;
yTest = yTest/sqrt(var(yTest)); % rescale the input to unit variance
T = length(yTest);

opts.verbose = 1; % view plots of the fitting process
  opts.minT = 5000;
  opts.numIts = 10;
[Var1,Lam1,om] = fit_probSTFT(yTest,D,opts); % trains filters to
                                              % match the spectrum
% clean spectrogram computed here
ZTest = probFB(yTest,Lam1,Var1,om,0); 
ATest = abs(ZTest').^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reconstruction (i.e. sampling missing data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = 1;
gapLim = [10,2000];
NGaps = length(gapPos);
gaps = ceil(linspace(gapLim(1),gapLim(2),L));

cnt = 0;

for l=[1]
  
  % display progress to user
  cnt= cnt+1;
  disp(['Progress ',num2str(cnt),'/',num2str(L)]);
  
  % signal with missing data
  yGap = yTest;
  ind = [];
  for ng=1:NGaps
    ind = [ind,gapPos(ng)+[-ceil(gaps(l)/2):+ceil(gaps(l)/2)]];
  end

  yGap(ind) = 0;
  
  ZGap = probFB(yGap,Lam1,Var1,om,0); 
  AGap = abs(ZGap').^2;
     
  varys = 1e-5*ones(T,1);
  varys(ind) = 1e5;
  
  % denoise using untrained (?) GTF
  verbose = 1;
  filter_only = 1;
%   tic
  ZRecon = probFB(yGap,Lam1,Var1,om,varys,verbose,filter_only);
%   toc
  yReconGTF = sum(real(ZRecon'),2);
  %ZTemp = probFB(yReconGTF,Lam1,Var1,om,0); 
  %AReconGTF = abs(ZTemp').^2;
  AReconGTF = abs(ZRecon').^2;
  
end

figure(2); clf
  subplot(3,1,1)
  plot(yTest)
  title('signal')

  subplot(3,1,2)
  plot(yReconGTF)
  title('reconstruction')

  subplot(3,1,3)
  plot(yGap)
  title('with gaps')
  
figure(3); clf
  subplot(3,1,1)
  imagesc(log(ATest)')
  title('signal')
  
  subplot(3,1,2)
  imagesc(log(AReconGTF)')
  title('reconstruction')
  
  subplot(3,1,3)
  imagesc(log(AGap)')
  title('with gaps')


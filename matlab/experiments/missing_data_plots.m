load('../../data/missing_data_results_speech.mat');
load('../../data/missing_data_results_speech_matern52.mat');
load('../../data/missing_data_results_speech_matern32_fast.mat');
load('../../data/missing_data_results_speech_matern52_fast.mat');

snr_y = missing_data_results.snr_y;
snr_y52 = missing_data_results_matern52.snr_y(:,:,3);
snr_y_fast = missing_data_results_matern32_fast.snr_y(:,:,2);
snr_y32_fast = missing_data_results_matern32_fast.snr_y(:,:,3);
snr_y52_fast = missing_data_results_matern52_fast.snr_y(:,:,3);
pesq_y = missing_data_results.pesq_y;
pesq_y52 = missing_data_results_matern52.pesq_y(:,:,3);
pesq_y_fast = missing_data_results_matern32_fast.pesq_y(:,:,2);
pesq_y32_fast = missing_data_results_matern32_fast.pesq_y(:,:,3);
pesq_y52_fast = missing_data_results_matern52_fast.pesq_y(:,:,3);
gaps = missing_data_results.gaps;
fs = missing_data_results.fs;
L = missing_data_results.L;
gapPos = missing_data_results.gapPos;
yTest = missing_data_results.yTest; % actual sig
yRecon1 = missing_data_results.yRecon1; % Matern 1/2
yRecon2 = missing_data_results.yRecon2; % Matern 3/2
yRecon3 = missing_data_results_matern52.yRecon2; % Matern 5/2

snr_y_all = zeros(10,L,7);
pesq_y_all = zeros(10,L,7);
snr_y_all(:,:,1:3) = snr_y;
snr_y_all(:,:,4) = snr_y52;
snr_y_all(:,:,5) = snr_y_fast;
snr_y_all(:,:,6) = snr_y32_fast;
snr_y_all(:,:,7) = snr_y52_fast;
pesq_y_all(:,:,1:3) = pesq_y;
pesq_y_all(:,:,4) = pesq_y52;
pesq_y_all(:,:,5) = pesq_y_fast;
pesq_y_all(:,:,6) = pesq_y32_fast;
pesq_y_all(:,:,7) = pesq_y52_fast;

% remove zeros from pesq_y stats
pesq_y_nonzero = [];
for i=1:size(pesq_y_all,1)
    if pesq_y_all(i,1,1) > 0
        pesq_y_nonzero = [pesq_y_nonzero; pesq_y_all(i,:,:)];
    end
end

% average across all audio recordings
median_pesq_y = median(pesq_y_nonzero,1);
median_snr_y = median(snr_y_all,1);
stderror_pesq = std(pesq_y_nonzero,1) / sqrt(size(pesq_y_all,1));
stderror_snr = std(snr_y_all,1) / sqrt(size(snr_y_all,1));

figure(1); clf

subplot(2,2,1)
hold on
f1=fill([gaps*1000/fs,fliplr(gaps*1000/fs)],[median_pesq_y(1,:,2)-stderror_pesq(1,:,2),fliplr(median_pesq_y(1,:,2)+stderror_pesq(1,:,2))],[1. 0.5 0.5]);
f2=fill([gaps*1000/fs,fliplr(gaps*1000/fs)],[median_pesq_y(1,:,3)-stderror_pesq(1,:,3),fliplr(median_pesq_y(1,:,3)+stderror_pesq(1,:,3))],[0.5 0.5 1.]);
f3=fill([gaps*1000/fs,fliplr(gaps*1000/fs)],[median_pesq_y(1,:,4)-stderror_pesq(1,:,4),fliplr(median_pesq_y(1,:,4)+stderror_pesq(1,:,4))],[0.5 1. 0.5]);
% f4=fill([gaps*1000/fs,fliplr(gaps*1000/fs)],[median_pesq_y(1,:,5)-stderror_pesq(1,:,5),fliplr(median_pesq_y(1,:,5)+stderror_pesq(1,:,5))],[1. 0.5 0.5]);
% f5=fill([gaps*1000/fs,fliplr(gaps*1000/fs)],[median_pesq_y(1,:,6)-stderror_pesq(1,:,6),fliplr(median_pesq_y(1,:,6)+stderror_pesq(1,:,6))],[0.5 0.5 1.]);
% f6=fill([gaps*1000/fs,fliplr(gaps*1000/fs)],[median_pesq_y(1,:,7)-stderror_pesq(1,:,7),fliplr(median_pesq_y(1,:,7)+stderror_pesq(1,:,7))],[0.5 1. 0.5]);
set(f1,'edgecolor',[1. 0.9 0.9]); set(f2,'edgecolor',[0.9 0.9 1.]); set(f3,'edgecolor',[0.9 1. 0.9]);
% set(f4,'edgecolor',[1. 0.9 0.9]); set(f5,'edgecolor',[0.9 0.9 1.]); set(f6,'edgecolor',[0.9 1. 0.9]);
alpha(f1,0.2); alpha(f2,0.2); alpha(f3,0.2);% alpha(f4,0.2); alpha(f5,0.2); alpha(f6,0.2);
% title('waveform')
plot(gaps*1000/fs,median_pesq_y(:,:,2),'r-')
plot(gaps*1000/fs,median_pesq_y(:,:,3),'b-')
plot(gaps*1000/fs,median_pesq_y(:,:,4),'g-')
% plot(gaps*1000/fs,median_pesq_y(:,:,5),'r--')
% plot(gaps*1000/fs,median_pesq_y(:,:,6),'b--')
% plot(gaps*1000/fs,median_pesq_y(:,:,7),'g--')
% legend('Matern 1/2','Matern 3/2')
xlabel('missing data duration (ms)', 'interpreter', 'Latex')
ylabel('PESQ', 'interpreter', 'Latex') % score')
% set(gca, 'YScale', 'log')
ylim([-Inf Inf])

subplot(2,2,2)
hold on
% f1=fill([gaps*1000/fs,fliplr(gaps*1000/fs)],[median_snr_y(1,:,2)-stderror_snr(1,:,2),fliplr(median_snr_y(1,:,2)+stderror_snr(1,:,2))],[1. 0.5 0.5]);
% f2=fill([gaps*1000/fs,fliplr(gaps*1000/fs)],[median_snr_y(1,:,3)-stderror_snr(1,:,3),fliplr(median_snr_y(1,:,3)+stderror_snr(1,:,3))],[0.5 0.5 1.]);
% f3=fill([gaps*1000/fs,fliplr(gaps*1000/fs)],[median_snr_y(1,:,4)-stderror_snr(1,:,4),fliplr(median_snr_y(1,:,4)+stderror_snr(1,:,4))],[0.5 1. 0.5]);
f1=patch([gaps*1000/fs,fliplr(gaps*1000/fs)],[median_snr_y(1,:,2)-stderror_snr(1,:,2),fliplr(median_snr_y(1,:,2)+stderror_snr(1,:,2))],[1. 0.5 0.5]);
f2=patch([gaps*1000/fs,fliplr(gaps*1000/fs)],[median_snr_y(1,:,3)-stderror_snr(1,:,3),fliplr(median_snr_y(1,:,3)+stderror_snr(1,:,3))],[0.5 0.5 1.]);
f3=patch([gaps*1000/fs,fliplr(gaps*1000/fs)],[median_snr_y(1,:,4)-stderror_snr(1,:,4),fliplr(median_snr_y(1,:,4)+stderror_snr(1,:,4))],[0.5 1. 0.5]);
% f4=fill([gaps*1000/fs,fliplr(gaps*1000/fs)],[median_snr_y(1,:,5)-stderror_snr(1,:,5),fliplr(median_snr_y(1,:,5)+stderror_snr(1,:,5))],[1. 0.5 0.5]);
% f5=fill([gaps*1000/fs,fliplr(gaps*1000/fs)],[median_snr_y(1,:,6)-stderror_snr(1,:,6),fliplr(median_snr_y(1,:,6)+stderror_snr(1,:,6))],[0.5 0.5 1.]);
% f6=fill([gaps*1000/fs,fliplr(gaps*1000/fs)],[median_snr_y(1,:,7)-stderror_snr(1,:,7),fliplr(median_snr_y(1,:,7)+stderror_snr(1,:,7))],[0.5 1. 0.5]);
% set(f1,'edgecolor',[1. 0.9 0.9]); set(f2,'edgecolor',[0.9 0.9 1.]); set(f3,'edgecolor',[0.9 1. 0.9]);
f1.EdgeAlpha=0; f2.EdgeAlpha=0; f3.EdgeAlpha=0;
% set(f4,'edgecolor',[1. 0.9 0.9]); set(f5,'edgecolor',[0.9 0.9 1.]); set(f6,'edgecolor',[0.9 1. 0.9]);
alpha(f1,0.2); alpha(f2,0.2); alpha(f3,0.2); %alpha(f4,0.2); alpha(f5,0.2); alpha(f6,0.2);
% p7=plot(gaps*1000/fs,median_snr_y(:,:,2),'k--');
p1=plot(gaps*1000/fs,median_snr_y(:,:,2),'r-');
p2=plot(gaps*1000/fs,median_snr_y(:,:,3),'b-');
p3=plot(gaps*1000/fs,median_snr_y(:,:,4),'g-');
% p4=plot(gaps*1000/fs,median_snr_y(:,:,5),'r--');
% p5=plot(gaps*1000/fs,median_snr_y(:,:,6),'b--');
% p6=plot(gaps*1000/fs,median_snr_y(:,:,7),'g--');
% plot(gaps*1000/fs,mean_snr_y(1,:,2)+stderror_snr(1,:,2))
% plot(gaps*1000/fs,mean_snr_y(1,:,2)-stderror_snr(1,:,2))
legend([p1,p2,p3],{'Matern 1/2','Matern 3/2','Matern 5/2'}, 'interpreter', 'Latex')
xlabel('missing data duration (ms)', 'interpreter', 'Latex')
ylabel('SNR (dB)', 'interpreter', 'Latex')
ylim([-Inf Inf])

gapnum = 2;
grey = [0.6 0.6 0.6];
% darkgrey = [0.25 0.25 0.25];
t = linspace(1,length(yTest),length(yTest))*1000/fs;
ind1gap = gapPos(gapnum)+[-ceil(gaps(L)/2):+ceil(gaps(L)/2)];
ind1 = gapPos(gapnum)+[-2.5*ceil(gaps(L)/2):+2.5*ceil(gaps(L)/2)];
t = t - t(ind1(1));
normaliser = max(yTest(ind1));
subplot(2,1,2)
plot(t(ind1),yTest(ind1)/normaliser,'Color',grey)
hold on
% plot(t(ind1gap),yTest(ind1gap)/normaliser, 'Color',darkgrey)
title('Signal Reconstruction', 'interpreter', 'Latex')
plot(t(ind1gap),yRecon1(ind1gap)/normaliser, 'r-')
plot(t(ind1gap),yRecon2(ind1gap)/normaliser, 'b-')
plot(t(ind1gap),yRecon3(ind1gap)/normaliser, 'g-')
xlabel('time (ms)', 'interpreter', 'Latex')
ylim([-1.25, 1.25])
legend('ground truth', 'interpreter', 'Latex')

%% TikZ / pgfplots

%   blue = [0 0 0.7];
  blue = [0.3 0.5 0.8];
%   red = [0.7 0 0];
  red = [0.9 0.4 0.2];
%   green = [0 0.7 0];
  green = [0.4 0.8 0.2];

close all
addpath ../matlab2tikz/src

figure(1); clf
hold on
f1=patch([gaps*1000/fs,fliplr(gaps*1000/fs)],[median_pesq_y(1,:,2)-stderror_pesq(1,:,2),fliplr(median_pesq_y(1,:,2)+stderror_pesq(1,:,2))],red);
f2=patch([gaps*1000/fs,fliplr(gaps*1000/fs)],[median_pesq_y(1,:,3)-stderror_pesq(1,:,3),fliplr(median_pesq_y(1,:,3)+stderror_pesq(1,:,3))],blue);
f3=patch([gaps*1000/fs,fliplr(gaps*1000/fs)],[median_pesq_y(1,:,4)-stderror_pesq(1,:,4),fliplr(median_pesq_y(1,:,4)+stderror_pesq(1,:,4))],green);
% f4=fill([gaps*1000/fs,fliplr(gaps*1000/fs)],[median_pesq_y(1,:,5)-stderror_pesq(1,:,5),fliplr(median_pesq_y(1,:,5)+stderror_pesq(1,:,5))],[1. 0.5 0.5]);
% f5=fill([gaps*1000/fs,fliplr(gaps*1000/fs)],[median_pesq_y(1,:,6)-stderror_pesq(1,:,6),fliplr(median_pesq_y(1,:,6)+stderror_pesq(1,:,6))],[0.5 0.5 1.]);
% f6=fill([gaps*1000/fs,fliplr(gaps*1000/fs)],[median_pesq_y(1,:,7)-stderror_pesq(1,:,7),fliplr(median_pesq_y(1,:,7)+stderror_pesq(1,:,7))],[0.5 1. 0.5]);
f1.EdgeAlpha=0; f2.EdgeAlpha=0; f3.EdgeAlpha=0;
% set(f4,'edgecolor',[1. 0.9 0.9]); set(f5,'edgecolor',[0.9 0.9 1.]); set(f6,'edgecolor',[0.9 1. 0.9]);
alpha(f1,0.15); alpha(f2,0.15); alpha(f3,0.15);% alpha(f4,0.2); alpha(f5,0.2); alpha(f6,0.2);
% title('waveform')
plot(gaps*1000/fs,median_pesq_y(:,:,2),'Color',red,'LineWidth',0.65)
plot(gaps*1000/fs,median_pesq_y(:,:,3),'Color',blue,'LineWidth',0.65)
plot(gaps*1000/fs,median_pesq_y(:,:,4),'Color',green,'LineWidth',0.65)
% plot(gaps*1000/fs,median_pesq_y(:,:,5),'r--')
% plot(gaps*1000/fs,median_pesq_y(:,:,6),'b--')
% plot(gaps*1000/fs,median_pesq_y(:,:,7),'g--')
% legend('Matern 1/2','Matern 3/2')
xlabel('missing data duration (ms)', 'interpreter', 'Latex', 'fontsize',14)
title('PESQ', 'interpreter', 'Latex','FontWeight','normal', 'fontsize',14) % score')
ylim([-Inf Inf])
box on

  % Save figure
  if true
    matlab2tikz('../../paper/fig/missing_data_matern12.tex', ...
      'noSize',true, ... 
      'relativeDataPath','./fig/', ...
      'extraAxisOptions',{'width=\figurewidth','height=\figureheight'},...
      'parseStrings',false, ...
      'checkForUpdates',false)  
  end


figure(2); clf
hold on
f1=patch([gaps*1000/fs,fliplr(gaps*1000/fs)],[median_snr_y(1,:,2)-stderror_snr(1,:,2),fliplr(median_snr_y(1,:,2)+stderror_snr(1,:,2))],red);
f2=patch([gaps*1000/fs,fliplr(gaps*1000/fs)],[median_snr_y(1,:,3)-stderror_snr(1,:,3),fliplr(median_snr_y(1,:,3)+stderror_snr(1,:,3))],blue);
f3=patch([gaps*1000/fs,fliplr(gaps*1000/fs)],[median_snr_y(1,:,4)-stderror_snr(1,:,4),fliplr(median_snr_y(1,:,4)+stderror_snr(1,:,4))],green);
% f4=fill([gaps*1000/fs,fliplr(gaps*1000/fs)],[median_snr_y(1,:,5)-stderror_snr(1,:,5),fliplr(median_snr_y(1,:,5)+stderror_snr(1,:,5))],[1. 0.5 0.5]);
% f5=fill([gaps*1000/fs,fliplr(gaps*1000/fs)],[median_snr_y(1,:,6)-stderror_snr(1,:,6),fliplr(median_snr_y(1,:,6)+stderror_snr(1,:,6))],[0.5 0.5 1.]);
% f6=fill([gaps*1000/fs,fliplr(gaps*1000/fs)],[median_snr_y(1,:,7)-stderror_snr(1,:,7),fliplr(median_snr_y(1,:,7)+stderror_snr(1,:,7))],[0.5 1. 0.5]);
f1.EdgeAlpha=0; f2.EdgeAlpha=0; f3.EdgeAlpha=0;
% set(f4,'edgecolor',[1. 0.9 0.9]); set(f5,'edgecolor',[0.9 0.9 1.]); set(f6,'edgecolor',[0.9 1. 0.9]);
alpha(f1,0.15); alpha(f2,0.15); alpha(f3,0.15); %alpha(f4,0.2); alpha(f5,0.2); alpha(f6,0.2);
% p7=plot(gaps*1000/fs,median_snr_y(:,:,2),'k--');
p1=plot(gaps*1000/fs,median_snr_y(:,:,2),'Color',red,'LineWidth',0.65);
p2=plot(gaps*1000/fs,median_snr_y(:,:,3),'Color',blue,'LineWidth',0.65);
p3=plot(gaps*1000/fs,median_snr_y(:,:,4),'Color',green,'LineWidth',0.65);
% p4=plot(gaps*1000/fs,median_snr_y(:,:,5),'r--');
% p5=plot(gaps*1000/fs,median_snr_y(:,:,6),'b--');
% p6=plot(gaps*1000/fs,median_snr_y(:,:,7),'g--');
% plot(gaps*1000/fs,mean_snr_y(1,:,2)+stderror_snr(1,:,2))
% plot(gaps*1000/fs,mean_snr_y(1,:,2)-stderror_snr(1,:,2))
% legend([p1,p2,p3],{'Mat\''ern-$\nicefrac{1}{2}$','Mat\''ern-$\nicefrac{3}{2}$','Mat\''ern-$\nicefrac{5}{2}$'}, 'interpreter', 'Latex')
legend([p1,p2,p3],{'$\nu=\nicefrac{1}{2}$','$\nu=\nicefrac{3}{2}$','$\nu=\nicefrac{5}{2}$'}, 'interpreter', 'Latex')
xlabel('missing data duration (ms)', 'interpreter', 'Latex', 'fontsize',14)
title('SNR (dB)', 'interpreter', 'Latex','FontWeight','normal', 'fontsize',14)
ylim([-Inf Inf])
box on


  % Save figure
  if true
    matlab2tikz('../../paper/fig/missing_data_matern32.tex', ...
      'noSize',true, ... 
      'relativeDataPath','./fig/', ...
      'extraAxisOptions',{'width=\figurewidth','height=\figureheight'},...
      'parseStrings',false, ...
      'checkForUpdates',false)  
  end



gapnum = 2;
grey = [0.6 0.6 0.6];
% darkgrey = [0.25 0.25 0.25];
t = linspace(1,length(yTest),length(yTest))*1000/fs;
ind1gap = gapPos(gapnum)+[-ceil(gaps(L)/2):+ceil(gaps(L)/2)];
ind1 = gapPos(gapnum)+[-2.*ceil(gaps(L)/2):+2.*ceil(gaps(L)/2)];
t = t - t(ind1(1));
normaliser = max(yTest(ind1));

figure(3); clf
plot(t(ind1),yTest(ind1)/normaliser,'Color',grey)
hold on
% plot(t(ind1gap),yTest(ind1gap)/normaliser, 'Color',darkgrey)
%title('Signal Reconstruction')
plot(t(ind1gap),yRecon1(ind1gap)/normaliser, 'Color',red,'LineWidth',0.535)
plot(t(ind1gap),yRecon2(ind1gap)/normaliser, 'Color',blue,'LineWidth',0.535)
plot(t(ind1gap),yRecon3(ind1gap)/normaliser, 'Color',green,'LineWidth',0.535)
xlabel('time (ms)', 'interpreter', 'Latex', 'fontsize',14)
% ylabel('Audio signal', 'interpreter', 'Latex')
title('Audio signal reconstruction', 'interpreter', 'Latex','FontWeight','normal', 'fontsize',14)
ylim([-1.25, 1.25])
%legend('ground truth')


  % Save figure
  if true
    matlab2tikz('../../paper/fig/missing_data.tex', ...
      'noSize',true, ... 
      'relativeDataPath','./fig/', ...
      'extraAxisOptions',{'width=\figurewidth','height=\figureheight'},...
      'parseStrings',false, ...
      'checkForUpdates',false)  
  end
clear; close all;
fs = 16000;
D = 40; % number of filters / frequency channels
kernel = 'matern32';
gapLim = [10, 320]; % 1-20ms
L = 7; % test this many different gap durations (use 7 for final results)
slow = 1; % 1 = Kalman smoother (slow & accurate), 0 = not implemented (coming soon...)
snr_y = zeros(10,L,3); pesq_y = zeros(10,L,3);
verbose = 1; % whether to show optimisation plots
rng('default');
rng(1,'twister');

%%
gapPos = [1500,5000,7000,9000,13000,18000]; % set manually to non-silent regions
gapPos_for_plot = gapPos;
[snr_y_exp, pesq_y_exp, yTest, yRecon1, yRecon2] = run_missing_data_experiment('../../audio/speech/speech0_female.wav',fs,D,kernel,L,gapLim,gapPos,slow,verbose);
snr_y(1,:,:)= snr_y_exp; pesq_y(1,:,:)= pesq_y_exp;

%%
gapPos = [500,1500,3500,5000,8000,10000]; % set manually to non-silent regions
[snr_y_exp, pesq_y_exp] = run_missing_data_experiment('../../audio/speech/speech1_male.wav',fs,D,kernel,L,gapLim,gapPos,slow,verbose);
snr_y(2,:,:)= snr_y_exp; pesq_y(2,:,:)= pesq_y_exp;

%%
gapPos = [1000,2500,4000,5500,6000,7500]; % set manually to non-silent regions
[snr_y_exp, pesq_y_exp] = run_missing_data_experiment('../../audio/speech/speech2_male.wav',fs,D,kernel,L,gapLim,gapPos,slow,verbose);
snr_y(3,:,:)= snr_y_exp; pesq_y(3,:,:)= pesq_y_exp;

%%
gapPos = [800,2200,5000,6500,10000,12500]; % set manually to non-silent regions
[snr_y_exp, pesq_y_exp] = run_missing_data_experiment('../../audio/speech/speech3_male.wav',fs,D,kernel,L,gapLim,gapPos,slow,verbose);
snr_y(4,:,:)= snr_y_exp; pesq_y(4,:,:)= pesq_y_exp;

%%
gapPos = [700,1600,2500,6000,8000,11000]; % set manually to non-silent regions
[snr_y_exp, pesq_y_exp] = run_missing_data_experiment('../../audio/speech/speech4_male.wav',fs,D,kernel,L,gapLim,gapPos,slow,verbose);
snr_y(5,:,:)= snr_y_exp; pesq_y(5,:,:)= pesq_y_exp;

%%
gapPos = [700,2000,3000,4000,5000,7000]; % set manually to non-silent regions
[snr_y_exp, pesq_y_exp] = run_missing_data_experiment('../../audio/speech/speech5_male.wav',fs,D,kernel,L,gapLim,gapPos,slow,verbose);
snr_y(6,:,:)= snr_y_exp; pesq_y(6,:,:)= pesq_y_exp;

%%
gapPos = [800,2000,3000,4000,10000,11000]; % set manually to non-silent regions
[snr_y_exp, pesq_y_exp] = run_missing_data_experiment('../../audio/speech/speech6_female.wav',fs,D,kernel,L,gapLim,gapPos,slow,verbose);
snr_y(7,:,:)= snr_y_exp; pesq_y(7,:,:)= pesq_y_exp;

%%
gapPos = [700,2000,5000,6000,9000,10000]; % set manually to non-silent regions
[snr_y_exp, pesq_y_exp] = run_missing_data_experiment('../../audio/speech/speech7_female.wav',fs,D,kernel,L,gapLim,gapPos,slow,verbose);
snr_y(8,:,:)= snr_y_exp; pesq_y(8,:,:)= pesq_y_exp;

%%
gapPos = [1000,2000,5000,8000,12000,13000]; % set manually to non-silent regions
[snr_y_exp, pesq_y_exp] = run_missing_data_experiment('../../audio/speech/speech8_female.wav',fs,D,kernel,L,gapLim,gapPos,slow,verbose);
snr_y(9,:,:)= snr_y_exp; pesq_y(9,:,:)= pesq_y_exp;

%%
gapPos = [1000,3500,7500,8500,10000,15000]; % set manually to non-silent regions
[snr_y_exp, pesq_y_exp] = run_missing_data_experiment('../../audio/speech/speech9_female.wav',fs,D,kernel,L,gapLim,gapPos,slow,verbose);
snr_y(10,:,:)= snr_y_exp; pesq_y(10,:,:)= pesq_y_exp;



%%
close all;
gaps = ceil(linspace(gapLim(1),gapLim(2),L));
% means across all audio recordings
mean_pesq_y = mean(pesq_y,1);
mean_snr_y = mean(snr_y,1);

figure(1); clf
subplot(2,1,1)
hold on
title('waveform')
plot(gaps*1000/fs,mean_pesq_y(:,:,2),'k--')
plot(gaps*1000/fs,mean_pesq_y(:,:,3),'k-')
legend('exp',kernel)
xlabel('gap /ms')
ylabel('PSEQ score')

subplot(2,1,2)
hold on
plot(gaps*1000/fs,mean_snr_y(:,:,2),'k--')
plot(gaps*1000/fs,mean_snr_y(:,:,3),'k-')
legend('exp',kernel)
xlabel('gap /ms')
ylabel('SNR')


figure(2); clf
subplot(2,1,1)
hold on
title('waveform')
plot(gaps*1000/fs,mean_pesq_y(:,:,2)-mean_pesq_y(:,:,1),'--k')
plot(gaps*1000/fs,mean_pesq_y(:,:,3)-mean_pesq_y(:,:,1),'-k')
legend('exp',kernel)
xlabel('gap /ms')
ylabel('PSEQ improvement')

subplot(2,1,2)
hold on
plot(gaps*1000/fs,mean_snr_y(:,:,2)-mean_snr_y(:,:,1),'--k')
plot(gaps*1000/fs,mean_snr_y(:,:,3)-mean_snr_y(:,:,1),'-k')
legend('exp',kernel)
xlabel('gap /ms')
ylabel('SNR improvement')


%%
missing_data_results_matern32 = struct;
missing_data_results_matern32.snr_y = snr_y;
missing_data_results_matern32.pesq_y = pesq_y;
missing_data_results_matern32.gaps = gaps;
missing_data_results_matern32.fs = fs;
missing_data_results_matern32.L = L;
missing_data_results_matern32.gapPos = gapPos_for_plot;
missing_data_results_matern32.yTest = yTest;
missing_data_results_matern32.yRecon1 = yRecon1;
missing_data_results_matern32.yRecon2 = yRecon2;
save('../../data/missing_data_results_speech_matern32.mat','missing_data_results_matern32');

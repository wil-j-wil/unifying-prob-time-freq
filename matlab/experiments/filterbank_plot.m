
%%
  fs = 2000;
  omega12 = pi ./ 3 - 0.35;
  omega32 = pi ./ 3;
  omega52 = pi ./ 3 + 0.35;
  
  % cts params
  var12 = 0.5;
  len = 5;  % must be > 1 otherwise lam_d <= 0
  lam12 = 1 ./ len;
  lam32 = 1.03 * sqrt(3) ./ len;
  lam52 = 1.03 * sqrt(5) ./ len;  
  
  % covariance and spectral density
  k12 = @(r) var12 * cos(omega12*r) * exp(-abs(r) * lam12);
  k32 = @(r) var12 * cos(omega32*r) * (1 + lam32*abs(r)) * exp(-abs(r) * lam32);
  k52 = @(r) var12 * cos(omega52*r) * (1 + lam52*abs(r) + (1/3)*lam52^2*r^2) * exp(-abs(r) * lam52);
  SD12 = @(w) var12 * lam12 * ((lam12^2 + (w - omega12)^2)^-1 + (lam12^2 + (w + omega12)^2)^-1);
  SD32 = @(w) 2 * var12 * lam32^3 * ((lam32^2 + (w - omega32)^2)^-2 + (lam32^2 + (w + omega32)^2)^-2);
  SD52 = @(w) (8/3) * var12 * lam52^5 * ((lam52^2 + (w - omega52)^2)^-3 + (lam52^2 + (w + omega52)^2)^-3);
  
  
  N = 1000;
  w_ = linspace(0,pi,N);
  r_ = linspace(-15,15,N);
  k12_vec = zeros(N,1); k32_vec = zeros(N,1); k52_vec = zeros(N,1);
  k12_mat = zeros(N,N); k32_mat = zeros(N,N); k52_mat = zeros(N,N);
  SD12_vec = zeros(N,1); SD32_vec = zeros(N,1); SD52_vec = zeros(N,1);
  for i=1:N
      k12_vec(i) = k12(r_(i));
      k32_vec(i) = k32(r_(i));
      k52_vec(i) = k52(r_(i));
      SD12_vec(i) = SD12(w_(i));
      SD32_vec(i) = SD32(w_(i));
      SD52_vec(i) = SD52(w_(i));
      for j=1:N
          k12_mat(i,j) = k12((r_(i)-r_(j))/3);
          k32_mat(i,j) = k32((r_(i)-r_(j))/3);
          k52_mat(i,j) = k52((r_(i)-r_(j))/3);
      end
  end
  k_mat = zeros(3*N,3*N);
  k_mat(1:N,1:N) = k12_mat;
  k_mat(N+1:2*N,N+1:2*N) = k32_mat;
  k_mat(2*N+1:end,2*N+1:end) = k52_mat;
  
  mu = zeros(N,1);
  rng('default');
  rng(7,'twister'); % 1, 2, 7 are good
  sample12 = mvnrnd(mu, k12_mat);
  sample32 = mvnrnd(mu, k32_mat);
  sample52 = mvnrnd(mu, k52_mat);
  
%   blue = [0 0 0.7];
  blue = [0.3 0.5 0.8];
%   red = [0.7 0 0];
  red = [0.9 0.4 0.2];
%   green = [0 0.7 0];
  green = [0.4 0.8 0.2];
  figure(1); clf
  subplot(2,2,1)
  hold on
  plot(w_*fs/(2*pi),SD12_vec,'Color',red,'LineWidth',1.75)
  plot(w_*fs/(2*pi),SD32_vec,'Color',blue,'LineWidth',1.75)
  plot(w_*fs/(2*pi),SD52_vec,'Color',green,'LineWidth',1.75)
  legend('Mat\''ern-$\nicefrac{1}{2}$','Mat\''ern-$\nicefrac{3}{2}$','Mat\''ern-$\nicefrac{5}{2}$')
  title('Filter response / Spectral density', 'interpreter', 'Latex')
  xlabel('Frequency (Hz)', 'interpreter', 'Latex')
  ylabel('Density (log scale)', 'interpreter', 'Latex')
  set(gca,'Yscale','log')
  ylim([10e-3, 8])
  
  subplot(2,2,2)
  hold on
  plot((r_)*1000/fs,k12_vec,'Color',red,'LineWidth',1.75)
  plot((r_)*1000/fs,k32_vec,'Color',blue,'LineWidth',1.75)
  plot((r_)*1000/fs,k52_vec,'Color',green,'LineWidth',1.75)
  xlim([-Inf Inf])
  title('Sinusoidal bases / Kernel functions', 'interpreter', 'Latex')
  xlabel('Time (ms)', 'interpreter', 'Latex')
%   legend('Matern 1/2','Matern 3/2','Matern 5/2')
  
  subplot(2,2,3)
  hold on
  im=imagesc(flipud(k_mat),[-0.25 0.5]);
  xlim([-Inf Inf])
  ylim([-Inf Inf])
%   colormap(fake_parula());
  colormap('default');
  axis('off')
  title('Covariance matrices', 'interpreter', 'Latex')
  im.AlphaData = .7;
  ax=gca;
%   caxis([-0.1 0.8]);
%   ax.CLim = [0 1];
  
  
  subplot(2,2,4)
  hold on
  plot((r_+15)*1000/fs,sample12,'Color',red,'LineWidth',1.2)
  plot((r_+15)*1000/fs,sample32,'Color',blue,'LineWidth',1.75)
  plot((r_+15)*1000/fs,sample52,'Color',green,'LineWidth',1.75)
  xlim([-Inf Inf])
  title('Freq. channel data / GP prior samples', 'interpreter', 'Latex')
  xlabel('Time (ms)', 'interpreter', 'Latex')
  
  %%
  
  figure(1); clf
  hold on
  plot(w_*fs/(2*pi),SD12_vec,'Color',red,'LineWidth',0.8)
  plot(w_*fs/(2*pi),SD32_vec,'Color',blue,'LineWidth',0.8)
  plot(w_*fs/(2*pi),SD52_vec,'Color',green,'LineWidth',0.8)
%   legend('Mat\''ern-$\nicefrac{1}{2}$','Mat\''ern-$\nicefrac{3}{2}$','Mat\''ern-$\nicefrac{5}{2}$')
  lg=legend('$\nu=\nicefrac{1}{2}$','$\nu=\nicefrac{3}{2}$','$\nu=\nicefrac{5}{2}$');
  
  
  title('Filter response / Spectral density', 'interpreter', 'Latex','FontWeight','normal')
  xlabel('frequency (Hz)', 'interpreter', 'Latex')
%   ylabel('Density (log scale)', 'interpreter', 'Latex')
  ylabel('\phantom{FOO}')
  set(gca,'Yscale','log')
  ylim([10e-3, 8])
  box on
  
  
  matlab2tikz('../../paper/fig/spec_density.tex', ...
      'noSize',true, ... 
      'relativeDataPath','./fig/', ...
      'extraAxisOptions',{'width=\figurewidth','height=\figureheight'},...
      'parseStrings',false, ...
      'checkForUpdates',false)  
  
  
  figure(2); clf
  hold on
  plot((r_)*1000/fs,k12_vec,'Color',red,'LineWidth',0.8)
  plot((r_)*1000/fs,k32_vec,'Color',blue,'LineWidth',0.8)
  plot((r_)*1000/fs,k52_vec,'Color',green,'LineWidth',0.8)
  xlim([-7.5 7.5])
  ylim([-0.5 0.52])
  title('Sinusoidal bases / Kernel functions', 'interpreter', 'Latex','FontWeight','normal')
  xlabel('time (ms)', 'interpreter', 'Latex')
%   legend('Matern 1/2','Matern 3/2','Matern 5/2')
  set(gca,'YTickLabel',[-1, 0, 1]);
  set(gca,'ytick',[-1, 0, 1]);
  box on

  matlab2tikz('../../paper/fig/kernels.tex', ...
      'noSize',true, ... 
      'relativeDataPath','./fig/', ...
      'extraAxisOptions',{'width=\figurewidth','height=\figureheight'},...
      'parseStrings',false, ...
      'checkForUpdates',false)  
  

  figure(3); clf
  hold on
  im=imagesc(flipud(k_mat),[-0.25 0.5]);
  xlim([-Inf Inf])
  ylim([-Inf Inf])
%   colormap(fake_parula());
  colormap('default');
  %axis('off')
  axis tight
  title('Covariance matrices', 'interpreter', 'Latex','FontWeight','normal')
  im.AlphaData = .7;
  ax=gca;
  box on
  set(gca,'XTickLabel','\phantom{3000}','YTickLabel','\phantom{3000}')
  xticks([500 1500 2500])
  xticklabels({'$\nu=\nicefrac{1}{2}$','$\nu=\nicefrac{3}{2}$','$\nu=\nicefrac{5}{2}$'})
  yticks([300 1500 2700])
  yticklabels({'$\nu=\nicefrac{5}{2}$','$\nu=\nicefrac{3}{2}$','$\nu=\nicefrac{1}{2}$'})
  xlabel('\phantom{FOO}')
  ylabel('\phantom{FOO}')
%   ax.YAxis.FontSize=3;
%   caxis([-0.1 0.8]);
%   ax.CLim = [0 1];

  matlab2tikz('../../paper/fig/covariances.tex', ...
      'noSize',true, ... 
      'relativeDataPath','./fig/', ...
      'extraAxisOptions',{'width=\figurewidth','height=\figureheight'},...
      'parseStrings',false, ...
      'checkForUpdates',false)  

  figure(4); clf
  hold on
  plot((r_+15)*1000/fs,sample12,'Color',red,'LineWidth',0.5)
  plot((r_+15)*1000/fs,sample32,'Color',blue,'LineWidth',0.8)
  plot((r_+15)*1000/fs,sample52,'Color',green,'LineWidth',0.8)
  xlim([0 15])
  title('Freq. channel data / Sample trajectories', 'interpreter', 'Latex','FontWeight','normal')
  xlabel('time (ms)', 'interpreter', 'Latex')
  ylim([-2, 2])
  box on

  matlab2tikz('../../paper/fig/samples.tex', ...
      'noSize',true, ... 
      'relativeDataPath','./fig/', ...
      'extraAxisOptions',{'width=\figurewidth','height=\figureheight'},...
      'parseStrings',false, ...
      'checkForUpdates',false)  

  
  
  
  

%%
 
 kernel12 = 'exp';
 kernel32 = 'matern32';
 kernel52 = 'matern52';
 kernelRBF = 'se';
 len = 10;
 lam12 = 1 / len;
 lam32 = sqrt(3) / len;
 lam52 = sqrt(5) / len;
 omega = pi/100;
 var12 = 1;
 var32 = 1.;
 var52 = 1.;
 fs=2*pi;
 N = 100;
 
 [A12,Q12,~,~,K12,~] = get_disc_model(lam12,var12,omega,D,kernel12,4);
 [A32,Q32,~,~,K32,~] = get_disc_model(lam32,var32,omega,D,kernel32,4);
 [A52,Q52,~,~,K52,~] = get_disc_model(lam52,var52,omega,D,kernel52,4);
  
  start = 0;
  z12 = zeros(2*K12,N); z32 = zeros(2*K32,N); z52 = zeros(2*K52,N);
  z12(1,1) = start; z32(1,1) = start; z32(3,1) = 0; z52(1,1) = start; z52(3,1) = 0; z52(5,1) = 0;
  Q32(1:2,1:2) = 0; Q32(3:end,3:end) = 0.000001; Q52(1:4,1:4) = 0; Q52(5:end,5:end) = 0.000001;
  for i=2:N
    z12(:,i) = A12*z12(:,i-1) + Q12*randn(2*K12,1);
    z32(:,i) = A32*z32(:,i-1) + Q32*randn(2*K32,1);
    z52(:,i) = A52*z52(:,i-1) + Q52*randn(2*K52,1);
  end
  
  figure(1);clf
  subplot(2,1,1)
  hold on
%   plot(z12(1,:))
  plot(z32(1,:))
  plot(z52(1,:))
  title('real part')
  subplot(2,1,2)
  hold on
%   plot(z12(1,:),z12(2,:))
  plot(z32(1,:),z32(2,:))
  plot(z52(1,:),z52(2,:))

 %%

  
  psi1 = 0.99;
  psi2 = 0.99999999;
  R = [cos(omega) -sin(omega); sin(omega) cos(omega)];
  A1 = psi1*R;
  A2 = psi2*R;
  rho1 = .1;
  fs=2*pi;
  N = 1.5*fs/omega;
  z12 = zeros(2,N);
  z_smooth = zeros(2,N);
  for i=2:N
    z12(:,i) = A1*z12(:,i-1) + rho1^2*randn(2,1);
    z_smooth(:,i) = A2*z_smooth(:,i-1) + rho1^2*[0; randn(1,1)];
  end
%   plotN = 100;
  figure(1);clf
  subplot(2,2,1)
  plot(z12(1,:))
  title('real part')
  subplot(2,2,3)
  plot(z12(1,:),z12(2,:))
  subplot(2,2,2)
  plot(z_smooth(1,:))
  title('real part')
  subplot(2,2,4)
  plot(z_smooth(1,:),z_smooth(2,:))
  figure(2);clf
  hold on
  plot(z12(1,:),z12(2,:))
  plot(z_smooth(1,:),z_smooth(2,:))
  
  
  
 %%
  %   blue = [0 0 0.7];
  blue = [0.3 0.5 0.8];
%   red = [0.7 0 0];
  red = [0.9 0.4 0.2];
%   green = [0 0.7 0];
  green = [0.4 0.8 0.2];
  
  
  N = 1000;
  omega = 2*pi/N;
  psi1 = 0.7;
  psi2 = 0.8;
  R = [cos(omega) -sin(omega); sin(omega) cos(omega)];
  A1 = psi1*R;
  A2 = psi2*R;
  rho1 = 2.;
  rho2 = 1.;
  z1 = zeros(2,N); z2 = zeros(2,N);
  c1 = zeros(N-1,1); s1 = zeros(N-1,1); c2 = zeros(N-1,1); s2 = zeros(N-1,1);
  for i=2:N
    z1(:,i) = A1*z1(:,i-1) + rho1^2*[0; randn(1,1)];
    z2(:,i) = A2*z2(:,i-1) + rho2^2*[0; randn(1,1)];
    c1(i-1) = cos(omega*i) + z1(1,i);
    s1(i-1) = sin(omega*i);
    c2(i-1) = cos(omega*i) + z2(1,i);
    s2(i-1) = sin(omega*i);
  end
%   plotN = 100;
  figure(2);clf
  hold on
  plot(c1,s1,'Color',blue,'LineWidth',3)
  plot(c2,s2,'Color',red,'LineWidth',3)
  xlim([-1.4 1.4])
  ylim([-1.4 1.4])
  

 
  
  
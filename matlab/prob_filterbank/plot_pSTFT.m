function figH = plot_pSTFT(varx,om,lamx,varargin)
  
  % function figH = plot_pSTFT(varx,om,lamx)
  % or 
  % function figH = plot_pSTFT(varx,om,lamx,FSamp)
  % or
  % function figH = plot_pSTFT(varx,om,lamx,FSamp,logSpec)
  
  % Plots the spectra of the probabilistic spectrogram
  
  % INPUTS
  % varx = innovations variance, size [D,1] (D = number of processes)
  % om = centre frequencies, size [D,1]
  % lamx = dynamical parameters, size [D,1] 
  %
  % Optional Inputs:
  % FSamp = sample rate  
  % logSpec = 1 if log-spectra to be plotted
  % 
  % OUTPUTS
  % figH = figure handle
  
    
  if nargin>3
    FSamp = varargin{1};
  else
    FSamp=1;
  end

  if nargin>4
    logSpec = varargin{2};
  else
    logSpec=0;
  end

  D = length(varx);
  NumPoints = 4000;
  
  freqs = linspace(0,0.5,NumPoints);
  spec = get_pSTFT_spec(freqs,lamx,varx,om);
  if nargin>5
    figH = figure(varargin{3}); clf
  else
      figH = figure
  end
  hold on;
  
  for d=1:D,
    plot(freqs*FSamp,spec(d,:))
  end
  
  plot(freqs*FSamp,sum(spec),'-m','linewidth',3)
 
  if logSpec==1
    set(gca,'yscale','log')
  end
 %set(gca,'xscale','log')
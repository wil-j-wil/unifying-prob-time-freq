  function figH = plot_fitAR2FB(pg,cosCF,cosDF,mVar,CFHist,DFHist, ...
				mVarHist,Objs,ObjCur,minVar,varargin)
  
  % function figH = plot_fitAR2FB(pg,cosCF,cosDF,mVar,CFHist,DFHist, ...
  %                               mVarHist,Obj,ObjCur,minVar,varargin)
  %
  % For plotting the parameter fit, objective and parameter
  % evolution for AR(2) fitting.

  
  [D,numLevels] = size(mVarHist);
  c2f = sum(~isnan(CFHist(1,:)))-1;
  
  cols = linspace(0,1,D)';
  cols = [cols,zeros(D,1),cols(end:-1:1)];

  position = [440 273 1060 800];
  clf;
  figH = gcf;
  set(gcf,'position',position);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Plot results
  
  left = 0.1;
  right = 0.01;
  top =0.01;
  bottom = 0.1;
  heightTop = 0.5;
  hspace = 0.1;
  vspace1 = 0.1;
  vspace2 = 0.02;
  
  widthTop = 1-left-right;
  widthBot = (1-left-right-hspace)/2;
  heightBot = (1-top-heightTop-bottom-vspace1-vspace2)/2;
  
  posax1 = [left,1-top-heightTop,widthTop,heightTop];
  posax21 = [left,bottom+heightBot+vspace2,widthBot,heightBot];
  posax22 = [left+widthBot+hspace,bottom+heightBot+vspace2,widthBot,heightBot];

  posax31 = [left,bottom,widthBot,heightBot];
  posax32 = [left+widthBot+hspace,bottom,widthBot,heightBot];

  ax21=axes('position',posax21);
  hold on;
  plot([1:length(Objs)],Objs,'-b')
  set(gca,'xticklabel','')
  ylabel('obj')

   if nargin>10
    likeHO = varargin{1};
    numLikeHO = sum(~isnan(likeHO));
    plot(linspace(1,length(Objs),numLikeHO),likeHO(1:numLikeHO),'-r','linewidth',2)
   legend('training','test - held out')
   title(['Held out like: ',num2str(likeHO(numLikeHO))])
   end

  
  ax31=axes('position',posax31);
  hold on;
  plot([length(Objs)-length(ObjCur)+1:length(Objs)],ObjCur,'-k')
  xlabel('iteration number')
  ylabel('cur obj')
    
  ax32=axes('position',posax32);
  hold on;
  ind = 1:c2f+1;
  for d=1:D
    plot(ind-1,mVarHist(d,ind),'color',cols(d,:));
  end
  plot([ind(1)-1,ind(end)-1],min(minVar)*[1,1],'-k')
  set(gca,'yscale','log')
  ylabel('mVar')
  xlabel('iteration number')
  
  ax22=axes('position',posax22);
  hold on;
  ind = 1:c2f+1;
  for d=1:D
    plot(ind-1,CFHist(d,ind),'color',cols(d,:));
    plot(ind-1,CFHist(d,ind)+DFHist(d,ind)/2,'color',cols(d,:));
    plot(ind-1,CFHist(d,ind)-DFHist(d,ind)/2,'color',cols(d,:));
  end
  set(gca,'ylim',[0,1/2]);
  set(gca,'xticklabel','')
  ylabel('CF/DF')

  
  ax1=axes('position',posax1);
  hold on
  NumFreqs = 500;
  RngFreqs = [0,1/2];
  for d=1:D
    [freqs,spec] = getSpecAR2CFDFmVar(cosCF(d),cosDF(d),mVar(d), ...
					       NumFreqs,RngFreqs);
    if d==1
      sum_spec = spec;
    else
      sum_spec = sum_spec+spec;
    end
    plot(freqs,spec,'-r','linewidth',2);
  end
  
  plot(freqs,sum_spec,'-b','linewidth',2)
  hplot = plot(linspace(0,1/2,length(pg)),pg,'-k','linewidth',2);
  uistack(hplot,'bottom');
  set(gca,'yscale','log')
  drawnow;

  function figH = plot_fit_probSTFT(pg,mVar,om,lamx,mVarHist,omHist, ...
				    lamHist,Objs,ObjCur,minVar,varargin)
  
  % function figH = plot_fitAR2FB(pg,cosCF,cosDF,mVar,CFHist,DFHist, ...
  %                               mVarHist,Obj,ObjCur,minVar,varargin)
  %
  % For plotting the parameter fit, objective and parameter
  % evolution for probabilistic spectrogram fitting
  
  
  [D,numLevels] = size(mVarHist);
  varx = mVar.*(1-lamx.^2); % conditional variance from marginal
  
  
  % NaNs indicate iterations which are yet to take place
  c2f = sum(~isnan(omHist(1,:)))-1;
  
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

  % only plot the data objective if there isn't held out data
  if nargin<11
    plot([1:length(Objs)],Objs,'-b')
  end
  
  set(gca,'xticklabel','')
  ylabel('obj')

   if nargin>10
    likeHO = varargin{1};
    numLikeHO = sum(~isnan(likeHO));
    plot(linspace(1,length(Objs),numLikeHO),likeHO(1:numLikeHO),'-r','linewidth',2)
   legend('test - held out')
   str1 = ['Held out like: ',num2str(likeHO(numLikeHO))];
   end

   if nargin>11
    likeUnReg = varargin{2};
    numLikeUnReg = sum(~isnan(likeUnReg));
    plot(linspace(1,length(Objs),numLikeUnReg),likeUnReg(1:numLikeUnReg),'-b','linewidth',2)
   legend('test - held out','training un-regularised')
      str1 = [str1,'   ','training like: ',num2str(likeUnReg(numLikeUnReg))];
   end
   
   if nargin>10
     title(str1)
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
    [fmax,df,varMa] = probSpec2freq(omHist(d,ind),lamHist(d,ind), ...
				    ones(1,length(ind)));

    patx = ind-1; patx = [patx,patx(end:-1:1)];
    paty = [fmax+df,fmax(end:-1:1)-df(end:-1:1)];
    plot(ind-1,fmax+df,'color',cols(d,:));
    plot(ind-1,fmax-df,'color',cols(d,:));
    %patch(patx,paty,cols(d,:),'edgecolor','none')
    plot(ind-1,fmax,'color',cols(d,:));

  end
  set(gca,'ylim',[0,1/2]);
  set(gca,'xticklabel','')
  ylabel('centre-freq and bandwidth')

  
  ax1=axes('position',posax1);
  hold on
  NumFreqs = 1000;%max([500,length(pg)];
    
  freqs = linspace(0,0.5,NumFreqs);
  spec = get_pSTFT_spec(freqs,lamx,varx,om);
  
  
  for d=1:D
    plot(freqs,spec(d,:),'-r','linewidth',2);
  end

  plot(freqs,sum(spec),'-b','linewidth',2)
  hplot = plot(linspace(0,1/2,length(pg)),pg,'-k','linewidth',2);
  uistack(hplot,'bottom');
  set(gca,'yscale','log')
  drawnow;

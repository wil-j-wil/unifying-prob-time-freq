function [varx,lamx,om,Info] = fit_probSTFT_SD(y,D,kernel,varargin)

% function [varx,lamx,om,Info] = fit_probSTFT(y,D,kernel,Opts)
%
% Will Wilkinson 20181026
% Fits STFT to the spectrum using the (continuous) spectral density of the
% GP kernel. Based on Richard Turner's original code.
%
% Fits a probabilstic STFT to a signal y.
% x_{1,t,d} = lam_d x_{1,t-1,d} +\eta_{1,t,d} \varx_d^{1/2}
% x_{2,t,d} = lam_d x_{2,t-1,d} +\eta_{2,t,d} \varx_d^{1/2}
% \eta_{i,t,d} ~ \Norm(0,1)
% y_t = real(\sum_{d} exp(i om_d)*(x_{1,t,d}+i x_{2,t,d}))
%
% Fits the parameters (lam_d, varx_d,om_d) such that the spectra are
% matched to the target spectrum specified in specTar.
%
% Specifically, we integrate out all of the hidden variables which
% renders the time-domain fitting problem into a spectrum matching
% problem. The log-likelihood being,
%
% like(\theta) = -1/2 \sum_t [ log(specModel_t(\theta)) +
%                 specTar_t/specModel_t(\theta) ]
%
% As spectra of y (specTar) are noisy (they have high variance) the
% likelihood has lots of local optima. In order to alleviate this
% effect we use a coarse to fine procedure based on using Welch's
% periodogram to estimate the spectrum of the signal y, rather than
% the standard fft method. This smooths the spectrum, reducing the
% variance at the cost of introducing some bias, but importantly it
% prevents local optima in the likelihood. We reduce the degree of
% smoothing through the optimisation (i.e. this is a coarse to fine
% method).
%
% i.e. we define
% specTar_{t,n} = a matrix containing the spectra of sections of y, then
% mn_specTar_t = 1/N * \sum_n specTar_{t,n}
% we covert this to a power spectral density estimate by dividing
% by the frequency resolution.
% We then fit the mean spectrum to mn_specTar.
%
% An alternative theoretical perspective is that we split the signal
% into a bunch of short chunks. We then model each as an independent
% draw from a probabilistic spectrogram. The likelihood is therefore a sum
% of the likelihood's derived above:
%
% like(\theta) = -1/2 \sum_{t,n} [ log(specModel_t(\theta)) +
%                 specTar_{t,n}/specModel_t(\theta) ]
%
% In order to do the fitting 

% The function parameterises the AR processes using the marginal
% variance of the processes rather than the conditional (marginal
% variance = conditional variance/(1-lam^2)) since the objective is
% more simply behaved in this space and the constraints on the
% parameters are simpler to enforce.
%
% See Chapter 5 of my thesis (Statistical Models for Natural Sounds
% by R.E.Turner) for more details about AR(2) Filter banks.
%
%
% INPUTS
% y = the signal, size [T,1]
% D = number of AR(2) processes to fit 
% OPTIONAL INPUTS:
% Opts = set of options
% 
% OUTPUTS
% varx = innovations variance, size [D,1] (D = number of processes)
% om = centre frequencies, size [D,1]
% lam = dynamical parameters, size [D,1] 
% Info = structure containing information about the estimation
%        process including
%        Objs = objectives at each scale concatenated into a vector
%        ins = number of completed iterations at each scale


% Rescale and zero-mean y (we rescale the variances at the end to
% compensate)
y = y(:) - mean(y(:));
varSig = var(y);
y = y/sqrt(varSig);

if ~isfield(varargin{1},'theta_init')
  % Initialise the STFT parameters
    % uniformly over log-frequency and make them fairly broad
    mVar = ones(D,1)/D;
    FLim = [1/40,0.35]; % limits for the centre frequency initialisation
    dfFrac = 1/20; % width of the processes wrt centre-frequency
%     fmax = logspace(log10(FLim(1)),log10(FLim(2)),D)';
    fmax = linspace(FLim(1),FLim(2),D)';
    [om,lamx] = freq2probSpec(fmax,fmax*dfFrac,1);

  cvar_d = mVar .* (1 - lamx.^2); % conditional variance
  % convert hypers to continuous form so we can calculate spectral density
  lam_c = -log(lamx);
  lam_max = 0.4;
  lamLim = ones(D,1)*[0,lam_max];
%   lam_c = min(-log(lamx), lam_max-1e-5);
  cvar_c = cvar_d ./ (1 - exp(-2 .* lam_c));
  mVar_c = cvar_c ./ (1 - lam_c.^2);
%   omLim = ones(D,1)*[0.0,pi]; % limits on the centre frequency
else
  theta_init = varargin{1}.theta_init;
  cvar_c = theta_init(1:D);
  if strcmp(kernel,'exp')
    lam_c = theta_init(D+1:2*D);
  elseif strcmp(kernel,'matern32')
    lam_c = theta_init(D+1:2*D) .* sqrt(3);
  elseif strcmp(kernel,'matern52')
    lam_c = theta_init(D+1:2*D) .* sqrt(5);
  end
  om = theta_init(2*D+1:end);
%   cvar_d = mVar .* (1 - lamx.^2); % conditional variance
  % convert hypers to continuous form so we can calculate spectral density
  if ~isfield(varargin{1},'bandwidth_lim')
      bandwidth_lim = 2;
  else
      bandwidth_lim = varargin{1}.bandwidth_lim;
  end
  lam_max = min(lam_c .* bandwidth_lim,1-1e-5);
%   lam_c = min(-log(lamx), lam_max-1e-5);
%   cvar_c = cvar_d ./ (1 - exp(-2 .* lam_c));
  mVar_c = max(cvar_c ./ (1 - lam_c.^2),1e-3);
  lamLim = [zeros(D,1),lam_max];
%   omLim = [om-0.1,om+0.1]; % limits on the centre frequency
  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constraints on the variables
minVar = max(mVar_c/400,1e-5); % minimum marginal variance - might want to taper
omLim = ones(D,1)*[0.0,pi]; % limits on the centre frequency

% if strcmp(kernel,'matern52')
%     keyboard
%   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COARSE TO FINE PROCESS

if ~isfield(varargin{1},'numLevels')
  numLevels = 40; % number of levels
else
  numLevels = varargin{1}.numLevels;
end

if ~isfield(varargin{1},'numIts')
  numIts = 10; % number of iterations per level
else
  numIts = varargin{1}.numIts;
end

if ~isfield(varargin{1},'minT')
  minT = 200;%200; % minimum segment size to compute spectrum of
else
  minT = varargin{1}.minT;
end

if ~isfield(varargin{1},'maxT')
  maxT = 1000; % maximum segment size to compute spectrum of
else
  maxT = varargin{1}.maxT;
end




% Observation noise tapered
if ~isfield(varargin{1},'vary_an')
  vary_an = logspace(log10(1e-6),log10(1e-10),numLevels); % annealing settings
else
  vary_an = varargin{1}.vary_an;
end

if ~isfield(varargin{1},'bet')
  bet = logspace(log10(100),0,numLevels); % how strongly to encourage variances to be pruned
else
  bet = logspace(log10(varargin{1}.bet),0,numLevels);
%  bet = varargin{1}.bet;
end


		 
T = length(y);

numFreq = floor(logspace(log10(min([minT,T])),log10(min([maxT,T])),numLevels));
ovLp = floor(numFreq/10);

Objs = [];

% keep the historical values of the variables as they evolve
omHist = repmat(NaN,[D,numLevels+1]);
lamHist = repmat(NaN,[D,numLevels+1]);
mVarHist = repmat(NaN,[D,numLevels+1]);

omHist(:,1) = om;
lamHist(:,1) = lam_c;
mVarHist(:,1) = mVar_c;

% Compute the Held Out likelihood if the user has requested it
if nargin>3
  if isfield(varargin{1},'yHO')
    compHOLike = 1;
    
    % compute the full periodogram for the held out data
    yHO = varargin{1}.yHO;
    yHO = yHO(:)/sqrt(varSig);
    
    THO = length(yHO);
    likeHO = repmat(NaN,[numLevels,1]);
	
    [pgHO,varpg] = welchMethod(yHO,THO,0);
    pgHO = pgHO/(1/2/THO);
    
    [pgUR,varpg] = welchMethod(y,THO,0);
  
    if mod(THO,2)==0
      % if even
      specHO = [pgHO;pgHO(end-1:-1:2)];
    else
      % if odd
      specHO = [pgHO;pgHO(end:-1:2)];
    end
  else
    compHOLike = 0;
  end
else
  compHOLike = 0;
end

 
% compute the full periodogram for the unregularised objective
likeUnReg = repmat(NaN,[numLevels,1]);

[pgUR,varpg] = welchMethod(y,T,0);
pgUR = pgUR/(1/2/T);
    
if mod(T,2)==0
  % if even
  specUR = [pgUR;pgUR(end-1:-1:2)];
else
  % if odd
  specUR = [pgUR;pgUR(end:-1:2)];
end
    

% convert hypers to continuous form so we can calculate spectral density
%   disc_cts_func = str2func(strcat('params_disc_to_cts_',kernel));
%   disc_cts_func = str2func(strcat('params_disc_to_cts_','exp'));
%   [om, mVar_c, lam_c] = disc_cts_func(om, mVar, lamx);
%   lam_c = -log(lamx);
%   var_c = mVar ./ (1 - exp(-2 .* lam_c));
% keyboard
for c2f=1:numLevels

  % Estmate spectrum of y using Welch's periodogram
  [pg,varpg] = welchMethod(y,numFreq(c2f),ovLp(c2f));
  pg = pg/(1/2/numFreq(c2f)); % convert to power spectral density estimate

  %pg = pg*numFreq;

  if mod(T,2)==0
    % if even
    specTar = [pg;pg(end-1:-1:2)];
  else
    % if odd
    specTar = [pg;pg(end:-1:2)];
  end  
  
%   minVar = var_c ./ 2.;
%   minVar(:) = 0.001;
  
  % Constraints on variables
 
  % current settings of the variables
  theta = [log(mVar_c-minVar);...
	   log(om-omLim(:,1))-log(omLim(:,2)-om);...
	   log(lam_c-lamLim(:,1))-log(lamLim(:,2)-lam_c)]; 
  % Base line noise to avoid divide by zero problems
  vary = max(specTar)*vary_an(c2f);
  
  Obj_func = strcat('get_Obj_pSTFT_',kernel);
  Obj_func_handle = str2func(Obj_func);
  
  
  % Fit the spectrum using probabilistc spectrogram
  [theta,ObjCur,inCur] = minimize(theta,Obj_func,numIts, ...
				  vary, specTar,minVar,omLim, ...
				  lamLim,bet(c2f)*numFreq(c2f)/numFreq(c2f(1)));


  if compHOLike==1
    % Compute HO likelihood if asked for
    [likeHO(c2f,1),dObjTemp] = Obj_func_handle(theta,0,specHO,minVar, ...
						   omLim,lamLim,0);    
  end
  
  % Compute the real objective on the full data too without regularisation
  [likeUnReg(c2f,1),dObjTemp] = Obj_func_handle(theta,0,specUR,minVar, ...
							 omLim,lamLim,0);
  
  
  % Collect information about the current iteration
  Objs = [Objs;ObjCur];
  ins(c2f) = inCur;
  
  % Compute the current settings of the parameters
  mVar_c= minVar+exp(theta(1:D));
  om = omLim(:,1)+(omLim(:,2)-omLim(:,1))./(1+exp(-theta(D+1:2*D)));
  lam_c = lamLim(:,1)+(lamLim(:,2)-lamLim(:,1))./(1+exp(-theta(2*D+1:3*D)));
  
  
  % convert back to discrete params.
%   cts_disc_func = str2func(strcat('params_cts_to_disc_',kernel));
%   cts_disc_func = str2func(strcat('params_cts_to_disc_','exp'));
%   [om, mVar, lamx] = cts_disc_func(om, mVar_c, lam_c);
%   mVar = var_c .* (1 - exp(-2 .* lam_c));
%   lamx = exp(-lam_c);

    
  % Store historical values for the processes
  omHist(:,c2f+1) = om;
  lamHist(:,c2f+1) = lam_c;
  mVarHist(:,c2f+1) = mVar_c;

  % plot if requested by user
  if nargin>3
    if isfield(varargin{1},'verbose')
      if varargin{1}.verbose==1
	% plot the fit
	if compHOLike==1
	  figH =  plot_fit_probSTFT_kern_cts(pg,mVar_c,om,lam_c,mVarHist,omHist,lamHist, ...
				    Objs,ObjCur,minVar,kernel,likeHO,likeUnReg);  
	else
	  figH =  plot_fit_probSTFT_kern_cts(pg,mVar_c,om,lam_c,mVarHist,omHist,lamHist, ...
				    Objs,ObjCur,minVar,kernel);  
	end
      end
    end
  end

% % re-initialised the pruned components of the model - don't reassign
% on the last iteration

if isfield(varargin{1},'reassign')&c2f<numLevels-1
  
  if ~isfield(varargin{1},'reassignThresh')
    reassignThresh = 8;
  else
    reassignThresh = varargin{1}.reassignThresh;
  end
  
  if varargin{1}.reassign==1
    for d=1:D    
      if mVar_c(d)/minVar(d)<reassignThresh

	disp('moving pruned process')
	  varx = mVar_c.*(1-lam_c.^2);
	  
	  % remove current process
	  lamxCur = lam_c; lamxCur(d) = [];
	  varxCur = varx; varxCur(d) = [];
	  omCur = om; omCur(d) = [];

	  N = length(specTar);
	  freqs = linspace(0,0.5,ceil(N/2));
	
      specModHandle = str2func(strcat('get_pSTFT_spec_cts_',kernel));
	  specMod = sum(specModHandle(freqs,lamxCur,varxCur,omCur));
	  
	  % reassigning based on differences in spectra
	  %dspec = specTar(1:ceil(N/2))-specMod';
	  
	  % reassigning based on differences in log-spectra
	  dspec = log(specTar(1:ceil(N/2)))-log(specMod');
	    
	  [val,pos] = max(dspec);
	  
	  mVar_c(d) = 1/20;
	  om(d) = 2*pi*freqs(pos);
	  lam_c(d) = 0.05;
	  varx = mVar_c.*(1-lam_c.^2);
	  
	  % specNew = sum(get_pSTFT_spec(freqs,lamx,varx,om));
	  %
	  % figure
	  % hold on
	  % plot(freqs,specTar(1:ceil(N/2)),'-k')
	  % plot(freqs,specMod,'-b')
	  % plot(freqs,specNew,'-r')
	  % set(gca,'yscale','log')
	  % keyboard
      end
    end
  end
end

% The following merging heuristic didn't perform well - I think that
% the shrinkage of the marginal variance should, to some extent, be
% sufficient to discourage two processes doing the same thing

% if isfield(varargin{1},'merge')&c2f<numLevels-1

%   if varargin{1}.merge==1
%     % find the centre frequencies and bandwidths
%     [fmax,df] = probSpec2freq(om,lamx,ones(D,1));
%     [val,ind] = sort(fmax);
%     fmax = fmax(ind)*numFreq(c2f);
%     df = df(ind)*numFreq(c2f);
    
    
%       for d=2:D    
% 	if fmax(d)-fmax(d-1)<1/2 &  abs(df(d)-df(d-1))<1/2
	  
% 	  disp('merging pruned process')
% 	  mVar(d-1) = mVar(d)+mVar(d-1);
% 	  om(d-1) = (om(d) + om(d-1))/2;
	  
% 	  % reassigning deleted process
	  
% 	  varx = mVar.*(1-lamx.^2);
	  
% 	  lamxCur = lamx; %lamxCur(d) = [];
% 	  varxCur = varx; %varxCur(d) = [];
% 	  omCur = om; %omCur(d) = [];
	  
% 	  N = length(specTar);
% 	  freqs = linspace(0,0.5,ceil(N/2));
	  
% 	  specMod = sum(get_pSTFT_spec(freqs,lamxCur,varxCur,omCur));
	  
% 	  % reassigning based on differences in spectra
% 	  %dspec = specTar(1:ceil(N/2))-specMod';
	  
% 	  % reassigning based on differences in log-spectra
% 	  dspec = log(specTar(1:ceil(N/2)))-log(specMod');
	  
% 	  [val,pos] = max(dspec);
	  
% 	  mVar(d) = 1/10;
% 	  om(d) = 2*pi*freqs(pos);
% 	  lamx(d) = 0.95;
% 	  varx = mVar.*(1-lamx.^2);
	  
% 	  % specNew = sum(get_pSTFT_spec(freqs,lamx,varx,om));
% 	  % %
% 	  %  figure
% 	  %  hold on
% 	  %  plot(freqs,specTar(1:ceil(N/2)),'-k')
% 	  %  plot(freqs,specMod,'-b')
% 	  %  plot(freqs,specNew,'-r')
% 	  %  legend('data','old model','new model')
% 	  %  set(gca,'yscale','log')
% 	  %  keyboard
% 	end
%       end
%   end
% end

end

Info.Objs = Objs;
Info.ins = ins;

% Rescale the variance parameters so their sum is equal to that of
% the input signal
lamx = lam_c;
mVar = mVar_c;
rescale = varSig/sum(mVar);
varx = mVar.*(1-lamx.^2);
varx = rescale*varx;

if compHOLike==1
  Info.likeHO=likeHO;
end

Info.likeUnReg=likeUnReg;

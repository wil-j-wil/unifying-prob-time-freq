function [Lam,Var,Info] = fitAR2FB(y,D,varargin)

% function [Lam,Var,Info] = fitAR2FB(y)
% Optional inputs:
% function [Lam,Var,Info] = fitAR2FB(y,Opts)
%
% Fits an AR(2) filter bank to a signal y.
%
% Each AR(2) process is:
% x_{d,t} = \lambda_1 x_{d,t-1} + \lambda_2 x_{d,t-2} + \sigma_d \epsilon
% \epsilon ~ Norm(0,1)
%
% And the signal is the sum 
% y_t = \sum_d x_{d,t}
%
% We integrate out all of the hidden variables which renders the
% time-domain fitting problem into a spectrum matching problem:
%
% like(\theta) = -1/2 \sum_t [ log(specAR2_t(\theta)) +
%                 specTar_t/specAR2_t(\theta) ]
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
% draw from an AR(2) filter bank. The likelihood is therefore a sum
% of the likelihood's derived above:
%
% like(\theta) = -1/2 \sum_{t,n} [ log(specAR2_t(\theta)) +
%                 specTar_{t,n}/specAR2_t(\theta) ]
%
% In order to do the fitting it is important NOT to use the AR(2)
% parameters directly (i.e. \lamda_{1:D,1:2} and \sigma_{1:D} as
% defined above). This is because some choice of these parameters
% result in non-stationary processes for which the power-spectrum is
% not defined and the analytical expression breaks down. The
% constraints on the parameters necessary to restrict the model to
% stationary processes are quite complicated and so we use a different
% parameterisation which makes things simple. In particular, we
% convert to a representation in terms of the cosine of the centre
% frequency of the processes:
% cosCF = \lambda_1 (\lambda_2-1)/(4 * \lambda_2))
% The cosine-bandwidth of the processes:
% cosDF = sqrt( -(1+\lamba_1^2+\lamda_2^2)/\lambda_2-4 cosCF^2-2)
% And the marginal variance of the processes:
% mVar = \sigma^2 (1-\lambd_2)/(1-\lamba_1.^2-\lambda_2.^2-\lambda_2-\lambda_1.^2.*\lambda_2+\lambda_2.^3);
% We can then place constraints on these variables and the
% optimisation problem is much better behaved.
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
% Lam = dynamical parameters, size [D,2] (D = number of AR2 processes)
% Var = innovations variance, size [D,1]
% Info = structure containing information about the estimation
%        process including
%        Objs = objectives at each scale concatenated into a vector
%        ins = number of completed iterations at each scale


% Rescale and zero-mean y (we rescale the variances at the end to
% compensate)
y = y(:) - mean(y(:));
varSig = var(y);
y = y/sqrt(varSig);

% Initialise the AR(2) process parameters
% uniformly over log-frequency and make them fairly broad
mVar = ones(D,1)/D;
FLim = [1/50,0.3];
CF = logspace(log10(FLim(1)),log10(FLim(2)),D)';
DF = CF/5;

% Convert from sample rate to cosine units - these are the
% coordinates used for the optimisation
[cosCF,cosDF] = CFDF2cosCFDF(CF,DF);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COARSE TO FINE PROCESS
numLevels = 20; % number of levels
numIts = 30; % number of iterations per level
minT = 200;%200; % minimum segment size to compute spectrum of
maxT = 3000; % maximum segment size to compute spectrum of

% Constraints on the variables
minVar = mVar/400; % minimum marginal variance - might want to taper
%cDFmin = 1/1000;%1/1000; % minimum of cosDF - the bandwidths in cosine
%                 % space - we might also want to taper this

cDFmin_an = logspace(log10(1e-3),log10(1e-5),numLevels);		 
cosDF(cosDF<cDFmin_an(1)) = cDFmin_an(1);

vary_an = logspace(log10(1e-6),log10(1e-10),numLevels);
		 
bet = 10; % how strongly to encourage variances to be pruned
		 
T = length(y);

numFreq = floor(logspace(log10(min([minT,T])),log10(min([maxT,T])),numLevels));
ovLp = floor(numFreq/10);

Objs = [];

CFHist = repmat(NaN,[D,numLevels+1]);
DFHist = repmat(NaN,[D,numLevels+1]);
mVarHist = repmat(NaN,[D,numLevels+1]);

CFHist(:,1) = CF;
DFHist(:,1) = DF;
mVarHist(:,1) = mVar;

if nargin>2
  if isfield(varargin{1},'yHO')
    compHOLike = 1;
    yHO = varargin{1}.yHO;
    yHO = yHO(:)/sqrt(varSig);
    
    THO = length(yHO);
    likeHO = repmat(NaN,[numLevels,1]);
    [pgHO,varpg] = welchMethod(yHO,THO,0);
    pgHO = pgHO/(1/2/THO);
  
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

  % Constraints on variables
  deltaDF = 1/2-abs(cosCF)/2+cosDF/4;%(1-abs(cosCF))/2;
  deltaCFmin = cosCF+1-deltaDF;
  deltaCFmax = 1-cosCF-deltaDF;
  
  cosCFLim =  [-deltaCFmin+cosCF,deltaCFmax+cosCF];

  cDFmin = cDFmin_an(c2f);  
  cosDFmin = cDFmin*ones(D,1);
  cosDFmax = max([2*deltaDF,cosDFmin*1.1]')';
  cosDFLim = [cosDFmin,cosDFmax];%[1/500,1];

  theta = [log(mVar-minVar);...
	   log(cosCF-cosCFLim(:,1))-log(cosCFLim(:,2)-cosCF);...
	   log(cosDF-cosDFLim(:,1))-log(cosDFLim(:,2)-cosDF)]; 

  % Base line noise to avoid divide by zero problems
  vary = max(specTar)*vary_an(c2f);
  
  % Fit the spectrum using AR(2) processes
  [theta,ObjCur,inCur] = minimize(theta,'getSpecAR2Obj',numIts, ...
				  vary,bet*numFreq(c2f)/numFreq(c2f(1)), ...
				  specTar,minVar,cosCFLim,cosDFLim);
  
  if compHOLike==1
    % Compute HO likelihood if asked for
    [likeHO(c2f,1),dObjTemp] = getSpecAR2Obj(theta,vary,0,specHO,minVar, ...
						   cosCFLim,cosDFLim);
  end
  
  % Collect information about the current iteration
  Objs = [Objs;ObjCur];
  ins(c2f) = inCur;
  
  % Compute the current settings of the parameters
  mVar = minVar+exp(theta(1:D));
  cosCF = cosCFLim(:,1)+(cosCFLim(:,2)-cosCFLim(:,1))./(1+exp(-theta(D+1:2*D)));
  cosDF = cosDFLim(:,1)+(cosDFLim(:,2)-cosDFLim(:,1))./(1+exp(- ...
						  theta(2*D+1:3*D)));
  % Compute the AR process parameters
  [Lam,Var] = cosFreq2AR2(cosCF,cosDF,mVar);
    
  % Store historical values for the processes
  [CF,DF] = cosCFDF2CFDF(cosCF,cosDF);
  CFHist(:,c2f+1) = real(CF);
  DFHist(:,c2f+1) = real(DF);
  mVarHist(:,c2f+1) = mVar;

  % plot if requested by user
  if nargin>2
    if isfield(varargin{1},'verbose')
      if varargin{1}.verbose==1
	% plot the fit
	if compHOLike==1
	figH = plot_fitAR2FB(pg,cosCF,cosDF,mVar,CFHist,DFHist,mVarHist, ...
			       Objs,ObjCur,minVar,likeHO);  
	else
	  figH = plot_fitAR2FB(pg,cosCF,cosDF,mVar,CFHist,DFHist,mVarHist, ...
			       Objs,ObjCur,minVar);  
	end
      end
    end
  end

  % re-initialised the pruned components of the model
  freq = linspace(0,1/2,length(pg));
       
  for d=1:D    
      if mVar(d)/minVar(d)<2
       specs = getSpecAR2cosFreq(cosCF,cosDF,mVar,freq);

       [val,loc]=max((log(pg+vary)-log(sum(specs)'+vary)));
      
       mVar(d) = 1/D;
       cosCFnew = cos(2*pi*freq(loc));
       if cosCFnew>0
	 cosCF(d) = min([1-3*cDFmin,cosCFnew]);
       else
	 cosCF(d) = max([-1+3*cDFmin,cosCFnew]);
       end
%       cosDF(d) = min([1/500,(1-abs(cosCF(d)))/2]);
       cosDF(d) = min([1.5*cDFmin,(1-abs(cosCF(d)))/2]);
     end
   end

  
  % for n=1:

  %     theta = [log(mVar-minVar);...
  % 	   log(cosCF-cosCFLim(:,1))-log(cosCFLim(:,2)-cosCF);...
  % 	   log(cosDF-cosDFLim(:,1))-log(cosDFLim(:,2)-cosDF)]; 
  %     [likeHO(c2f,1),dObjTemp] = getSpecAR2Obj(theta,vary,0,specTar,minVar,cosCFLim,cosDFLim);
  % end

end

Info.Objs = Objs;
Info.ins = ins;

% Rescale the variance parameters so their sum is equal to that of
% the input signal
[CF,dF, mVar] = AR22freq(Lam,Var);
rescale = varSig/sum(mVar);
Var = rescale*Var;


if compHOLike==1
  Info.likeHO=likeHO;
end

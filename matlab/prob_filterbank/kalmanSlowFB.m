function  [lik,Xfin,Pfin,varargout] = kalmanSlowFB(lamx,varx,om, ...
						  vary,y,varargin)


% function
% [lik,Xfin,Pfin]=kalmanSlowFB(lamx,varx,om,vary,y);
% 
% Implements the probabilistic filter bank (see Cemgil and Godsill,
% 2005). Optionally returns the sufficient statistics of the Gaussian
% LDS. Based on Zoubin's code. Modified by Richard Turner.
%
% x_{t}|x_{t-1} ~ Norm(A x_{t-1},Q)
% y_{t}|x_{t} ~ Norm(C x_{t},R) 
% x_1 ~ Norm(x0,P0)  
%
% With optional outputs and inputs:
%
% function
% [lik,Xfin,Pfin,Ptsum,YX,A1,A2,A3]=kalmanSlowFB(lamx,varx,om,vary,Y,verbose,KF);
%
% see test_kalmanSlowFB.m for unit tests.
% 
% INPUTS:
% lamx = dynamical AR parameters [D,1]
% varx = dynamical noise parameters [D,1]
% om = mean frequencies of the sinusoids [D,1]
% vary = observation noise (either scalar or vector [T,1])
% y = Data, size [T,1] 
%
% OPTIONAL INPUTS:
% verbose = binary scalar, if set to 1 displays progress
%           information
% KF = binary scalar, if set to 1 carries out Kalman Filtering
%      rather than Kalman smoothing. Cannot return the sufficient
%      statistics in this case i.e. Ptsum, YX, A1, A2 and A3.
%
% OUTPUTS
% lik = likelihood
% Xfin = posterior means, size [N,K,T]
% Pfin = posterior covariance, size [K,K,T]
%
% OPTIONAL OUTPUTS:
% Ptsum = \sum_{t=1}^T <x_{k t} x_{k' t}>, size [K,K]
% YX = \sum_{t=1}^T y_{t}<x_{k t}>, size [D,K]
% A1 = \sum_{t=2}^T <x_{k t} x_{k' t-1}>, size [K,K] 
% A2 = \sum_{t=2}^T <x_{k t-1} x_{k' t-1}>, size [K,K]
% A3 = \sum_{t=2}^T <x_{k t} x_{k' t}>, size [K,K]
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K = length(lamx);
T = length(y);

C = repmat([1,0],[1,K]);

A = zeros(2*K);

for k=1:K
  ind = [1:2]+(k-1)*2;
  A(ind,ind) = lamx(k)*[cos(om(k)),-sin(om(k));sin(om(k)),cos(om(k))];
end

temp2 = [varx';
         varx'];

Q = diag(temp2(:));

T_R = length(vary);

x0 = zeros(2*K,1);

mVar = varx./(1-lamx.^2);
temp3 = [mVar';
	 mVar'];

P0 = diag(temp3(:));

Y = reshape(y,[1,1,T]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[N, D, T]=size(Y);
K=length(x0);
tiny=exp(-700);
I=eye(K);
const=(2*pi)^(-D/2);
problem=0;
lik=0;

Xpre=zeros(N,K);   % P(x_t | y_1 ... y_{t-1})
Xcur=zeros(N,K,T);   % P(x_t | y_1 ... y_t)
Xfin=zeros(N,K,T);   % P(x_t | y_1 ... y_T)    given all outputs

Ppre=zeros(K,K,T);
Pcur=zeros(K,K,T);
Pfin=zeros(K,K,T); 

Pt=zeros(K,K); 
Pcov=zeros(K,K); 
Kcur=zeros(K,D);
invP=zeros(D,D);
J=zeros(K,K,T);

if nargin<=5
  verbose = 0 ;
else
  verbose = varargin{1};
end

if nargin<=6
  KF = 0 ;
else
  KF = varargin{2};
end

if verbose==1
  if KF==1
    disp('Kalman Filtering')
  else
    disp('Kalman Smoothing')
  end
end
tic
%%%%%%%%%%%%%%%
% FORWARD PASS

%R=R+(R==0)*tiny;

Xpre=ones(N,1)*x0';
Ppre(:,:,1)=P0;

CntInt=T/5; % T / 2*number of progress values displayed

for t=1:T
  
  if T_R>1
    R = vary(t);
    invR = inv(R);
  else
    R = vary;
    invR = inv(R);
  end
  
  if verbose==1 && mod(t-1,CntInt)==0
    fprintf(['Progress ',num2str(floor(50*t/T)),'%%','\r'])
  end
    
  if (K<D)
    temp1= R\C;%  rdiv(C,R);
    temp2=temp1*Ppre(:,:,t); % inv(R)*C*Ppre
    temp3=C'*temp2;
   % temp4=inv(I+temp3)*temp1';
    temp4=(I+temp3)\temp1';

    invP=invR-temp2*temp4; 
    CP= temp1' - temp3*temp4;    % C'*invP
  else
%    temp1=diag(R)+C*Ppre(:,:,t)*C';
    temp1=R+C*Ppre(:,:,t)*C';
    invP=inv(temp1);
%    CP=C'*invP;
    CP=C'/temp1;
  end

  Kcur=Ppre(:,:,t)*CP;
  KC=Kcur*C;
  Ydiff=Y(:,:,t)-Xpre*C';
  Xcur(:,:,t)=Xpre+Ydiff*Kcur'; 
  %size(Xcur)
  Pcur(:,:,t)=Ppre(:,:,t)-KC*Ppre(:,:,t);
% keyboard
  if (t<T)
    Xpre=Xcur(:,:,t)*A';
    Ppre(:,:,t+1)=A*Pcur(:,:,t)*A'+Q;
  end

  % calculate likelihood
  
  % Old version of the code
  %  detiP=sqrt(det(invP));
  % if (isreal(detiP) & detiP>0)
  %   lik=lik+N*log(detiP)-0.5*sum(sum(Ydiff.*(Ydiff*invP)));
  % else
  %   problem=1;
  % end;

  logdetiP=sum(log(diag(chol(invP))));
  lik=lik+N*logdetiP-0.5*sum(sum(Ydiff.*(Ydiff*invP)));
  
end

lik=lik+N*T*log(const);

% Figure out whether the user wants to do Kalman Smoothing or
% Kalman Filtering
toc
if KF==1
  % Only Kalman filtering

  Xfin=Xcur;
  Pfin=Pcur; 
else
  % Kalman Filtering
  
  %%%%%%%%%%%%%%%
  % BACKWARD PASS
  
  t=T; 
  Xfin(:,:,t)=Xcur(:,:,t);
  Pfin(:,:,t)=Pcur(:,:,t); 
      
  for t=(T-1):-1:1
    if verbose==1 && mod(t+1,CntInt)==0
      fprintf(['Progress ',num2str(50+floor(50*(T-t+1)/T)),'%%','\r'])
    end
    
    J(:,:,t)=Pcur(:,:,t)*A'/Ppre(:,:,t+1);
    Xfin(:,:,t)=Xcur(:,:,t)+(Xfin(:,:,t+1)-Xcur(:,:,t)*A')*J(:,:,t)';
    Pfin(:,:,t)=Pcur(:,:,t)+J(:,:,t)*(Pfin(:,:,t+1)-Ppre(:,:,t+1))*J(:,:,t)';
    
  end
end

% If the user requests the additional sufficient statistics
if nargout>3 && KF==1

  disp('Error: Have not implemented sufficient statistics computation for Kalman Filtering')
  return;
else
  
  A1=zeros(K);
  A2=zeros(K);
  A3=zeros(K);
  Ptsum=zeros(K);
  YX=zeros(D,K);    
    
  Pt = Pfin(:,:,T) + Xfin(:,:,T)'*Xfin(:,:,T)/N; 
  A2 = -Pt;
  Ptsum = Pt;    
  YX = Y(:,:,T)'*Xfin(:,:,T);

  for t=(T-1):-1:1
    Pt=Pfin(:,:,t) + Xfin(:,:,t)'*Xfin(:,:,t)/N; 
    Ptsum=Ptsum+Pt;
    YX=YX+Y(:,:,t)'*Xfin(:,:,t);
  end
          
  A3 = Ptsum-Pt;
  A2 = Ptsum+A2;
    
  Pcov=(I-KC)*A*Pcur(:,:,T-1);
  A1=A1+Pcov+Xfin(:,:,T)'*Xfin(:,:,T-1)/N;
    
  for t=(T-1):-1:2      
    Pcov=(Pcur(:,:,t)+J(:,:,t)*(Pcov-A*Pcur(:,:,t)))*J(:,:,t-1)';
    A1=A1+Pcov+Xfin(:,:,t)'*Xfin(:,:,t-1)/N;
  end
    
  if problem 
    fprintf('!!!!!!!!! problem  !!!!!!!!!'); problem=0;
    lik = NaN;
  end
    
  % optional output arguments
  varargout(1) = {Ptsum};
  varargout(2) = {YX};
  varargout(3) = {A1};
  varargout(4) = {A2};
  varargout(5) = {A3};
end

if verbose==1
  fprintf('                                        \r')
end
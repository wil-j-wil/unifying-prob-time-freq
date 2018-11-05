function  [lik,Xfin,Pfin,varargout] = kernel_ss_kalmanSlowFB_rewrite(A,Q,C,P0,K,vary,y,varargin)
% 
%                      
%                     
%% Some of the tweaks from the original function
                      
T = length(y);
Y = reshape(y,[1,1,T]);
lik=0;
if length(vary) == 1
    vary = vary * ones(T,1);
end

if nargin<=7
  verbose = 0 ;
else
  verbose = varargin{1};
end

if nargin<=8
  KF = 0 ;
else
  KF = varargin{2};
end

if verbose==1
  if KF==1
    disp('Kalman Filtering (Arno''s re-written filter)')
  else
    disp('Kalman Smoothing (Arno''s re-written filter)')
  end
end

%% ARNO's re-write starts here

    tic

    % T / 2*number of progress values displayed
    CntInt=T/5; 

    % Set initial state and assign model matrices
    m = zeros(size(A,1),1);
    P = P0;
    H = C;
    
    % Allocate space for results
    MS = zeros(size(m,1),size(Y,3));
    PS = zeros(size(m,1),size(m,1),size(Y,3));
    
    % ### Forward filter
    
    disp('filtering stage:')
    
    % The filter recursion
    for k=1:size(Y,3)
        
        % Progress
        if verbose==1 && mod(k-1,CntInt)==0
            fprintf(['Progress ',num2str(floor(50*k/size(Y,3))),'%%','\r'])
        end
        
        R = vary(k);
        
        % Prediction
        if (k>1)
            m = A*m;
            P = A*P*A' + Q;
        end
        
        % Kalman update
        S = H*P*H' + R;
        K = P*H'/S;
        v = Y(:,:,k)-H*m;
        m = m + K*v;
        P = P - K*H*P;
      
        % Evaluate the energy (neg. log lik): Check this
        lik = lik + .5*size(S,1)*log(2*pi) + .5*log(S) + .5*v'/S*v;

        % Store estimate
        MS(:,k)   = m;
        PS(:,:,k) = P;
        
    end
  
    toc
       
    % ### Backward smoother

    % Should we run the smoother?
    if KF==1
        
        % Only return filter  
        
    else
        
        disp('smoothing stage:')
      
        % Rauch-Tung-Striebel smoother
        for k=size(MS,2)-1:-1:1
            
            % Progress
            if verbose==1 && mod(k+1,CntInt)==0
                fprintf(['Progress ',num2str(50+floor(50*(size(MS,2)-k+1)/size(MS,2))),'%%','\r'])
            end
            
            % Smoothing step (using Cholesky for stability)
            PSk = PS(:,:,k);
            
            % Pseudo-prediction
            PSkp = A*PSk*A'+Q;
            
            % Solve the Cholesky factorization
            [L,notposdef] = chol(PSkp,'lower');
            
            % Numerical problems in Cholesky, retry with jitter
            if notposdef>0
                jitterSigma2 = 1e-9;
                jitter = sqrt(jitterSigma2)*diag(rand(size(A,1),1));
                L = chol(PSkp+jitter,'lower');
            end
            
            % Continue smoothing step
            G = PSk*A'/L'/L;
            
            % Do update
            m = MS(:,k) + G*(m-A*MS(:,k));
            P = PSk + G*(P-PSkp)*G';
            
            % Store estimate
            MS(:,k)   = m;
            PS(:,:,k) = P;
            
        end

        toc
    end
    
    % Return variables of interest
    lik = -lik;
    Xfin = reshape(MS,[1 size(MS)]);
    Pfin = PS;
   
%%
    
   if verbose==1
     fprintf('                                        \r')
   end    
    

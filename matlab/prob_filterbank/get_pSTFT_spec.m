function spec = get_pSTFT_spec(freqs,lamx,varx,om)

% function spec = get_pSTFT_spec(freqs,lamx,varx,om)
%
% Computes the spectrum of a probabilistic STFT model
% x_{1,t} = lam x_{1,t-1} +\eta_{1,t} \varx^{1/2}
% x_{2,t} = lam x_{2,t-1} +\eta_{2,t} \varx^{1/2}
% \eta_{i,t} ~ \Norm(0,1)
% y_t = Real(exp(i om)*(x_{1,t}+i x_{2,t}))
%
% see getCompSpecPFB.m for the complex spectrum
%
% INPUTS
% freqs = frequencies at which to evaluate the spectrum, [1,N]
% lamx = dynamical AR parameters [D,1]
% varx = dynamical noise parameters [D,1]
% om = mean frequencies of the sinusoids [D,1]
%
% OUTPUTS
% spec = spectrum of the probabilistic STFT model
% 

D = length(lamx);
N = length(freqs);

spec = zeros(D,N);
omegas = 2*pi*freqs(:)';

for d=1:D

  spec(d,:) = 1/2*varx(d)./(1+lamx(d).^2-2*lamx(d)*cos(omegas-om(d))) ...
            + 1/2*varx(d)./(1+lamx(d).^2-2*lamx(d)*cos(omegas+om(d))); 
end



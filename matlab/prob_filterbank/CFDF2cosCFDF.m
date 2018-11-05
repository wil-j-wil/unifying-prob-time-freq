function [cosCF,cosDF] = CFDF2cosCFDF(CF,DF)

% function [cosCF,cosDF] = CFDF2cosCFDF(CF,DF)
%
% Converts centre frequency and bandwidth from frequency space (in
% units of samples) to cosine units.
%
% INPUTS
% CF = centre frequency in samples size [D,1]
% DF = bandwidth in samples size [D,1]
%
% OUTPUTS
% cosCF = cos centre frequencies size [D,1]
% cosDF = cos bandwidths size [D,1]

D = length(CF);

cosCF = cos(2*pi*CF);

%deltas = linspace(0,1/2,10000);
max_delta = 2*min([1-cosCF,1+cosCF]')';

cosDF = zeros(D,1);

for d=1:D
  deltas = linspace(0,max_delta(d),10000);
  obj=abs(DF(d)-(acos(cosCF(d)-deltas/2)-acos(cosCF(d)+deltas/2))/(2*pi));
%  plot(deltas,obj)
%  keyboard
  [val,pos] = min(obj);
  cosDF(d) = deltas(pos);
end


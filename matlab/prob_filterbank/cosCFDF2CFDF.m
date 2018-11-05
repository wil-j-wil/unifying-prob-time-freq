function [CF,DF] = cosCFDF2CFDF(cosCF,cosDF)

% function [cosCF,cosDF] = cosCFDF2CFDF(cosCF,cosDF)
%
% Converts cosine centre frequency and bandwidth into frequency space
% (in units of samples).
%
% INPUTS
% cosCF = cos centre frequencies size [D,1]
% cosDF = cos bandwidths size [D,1]
%
% OUTPUTS
% CF = centre frequency in samples size [D,1]
% DF = bandwidth in samples size [D,1]

D = length(cosCF);

CF = acos(cosCF)/(2*pi);

DF = 1/(2*pi)*(acos(cosCF-cosDF/2)-acos(cosCF+cosDF/2));

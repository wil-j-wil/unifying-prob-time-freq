function [x,frac] = loud_floor(x,xfloor)

% function x = loud_floor(x,xfloor)
%
% apply a floor to a signal

ind = x<xfloor;
frac = sum(ind(:))/length(x(:));
x(ind) = xfloor;


function [cc] = CC(n,niuu,x,alfa)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%cc=gamma(n*alfa+1)/gamma(alfa*(n+1))*x.^( alfa*(n+1)-1);

%cc=(-1)*(n*alfa).^2/gamma(niuu+2)*x.^(niuu+1)*hypergeom(1,[niuu/2+1 (1+niuu)/2+1],(-1)*(n*alfa*x).^2/4);
cc=gamma(n*alfa+1)/gamma(alfa*(n+1))*x.^( alfa*(n+1)-1);



end


function [ddd] = Dkrn(kk1,nn1,niuu,ii,alfarr)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%ddd=gamma(nn1*alfarr+1)/gamma(nn1*alfarr+niuu)*(ii-kk1).^(nn1*alfarr+niuu-1);

%ddd=(-1)*(nn1*alfa).^2/gamma(niuu+2)*(ii-nn1*alfa*kk1).^(niuu+1)*U(ii-nn1*alfa*kk1);
%ddd=ddd*hypergeom(1,[niuu/2+1 (1+niuu)/2+1],(-1)*(nn1*alfa).^2*(ii-nn1*alfa*kk1).^2/4);
ddd=gamma(nn1*alfarr+1)/gamma(nn1*alfarr+niuu)*(ii-kk1).^(nn1*alfarr+niuu-1);


end


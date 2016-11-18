function [ H ] = Heaviside( x , eps)
% Numerical (smooth) Heaviside function:
%                  { 1                                               x > eps
% H_eps(x) = { 0                                               x < - eps
%                  { 0.5*(1+(2/pi)*arctan(x/eps))       |x| <= eps

if nargin == 1
    eps = 1;
end

H = zeros(size(x));

H( abs(x) <= eps) = 0.5*(1+(2/pi)*atan(x(abs(x)<=eps)/eps)) ;
H( x > eps) = 1;
% H(x<-eps) = 0 by default

end


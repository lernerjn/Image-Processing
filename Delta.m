function [ D ] = Delta( x , eps)
% Numerical (smooth) Heaviside function:
%
% Delta_eps(x) = (1/pi) * ( eps/(eps^2+x^2))

if nargin == 1
    eps = 2;
end

D = (abs(x)<eps).*(1/pi) .* ( eps./ (eps^2 + x.^2));

end


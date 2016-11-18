% Kernel for the Bhattacharyaa Segmentation
% Note that we divide by normalization to allow the user to ensure that integral( K(z) dz) = 1
%
% Guassian Kernel:
%       K(y) = (1/(normalization*sqrt(2*pi*sigma^2))) * exp( -y^2 / ( 2*sigma^2))
%       where sigma^2 is the variance.
%
% Delta Function:
%       K(y) =  [ (1/(normalization*2*sigma))*(1+cos(pi*y/sigma))      |y| < sigma 
%                  [ 0                                                      otherwise

function [out] = Kernel(y, sigma, normalization)
    if nargin == 1
        sigma = 1;
        normalization =1;
    elseif nargin == 2
        normalization = 1;
    end
    % Gaussian Kernel (see between equations (2) and (3) in Seeing the Unseen)
    out =(1/(normalization*sqrt(2*pi*sigma^2))).* exp(  -(y).^2 ./(2*sigma^2))  ;
    % Delta Function see equation (6) in Seeing the Unseen
%     out = ( abs(y) < sigma) .* 1/(normalization*2*sigma) .* ( 1 + cos(pi*y/sigma));
end

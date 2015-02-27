function [ I ] = HeatEquationConv( I, tend, p, strIn)
%Author: Jeremy Lerner, Stony Brook University
%  [ Io ] = HeatEquationConv( J, 0, 1, 1, 'Conv Ben' );

%This program runs the non linear heat equation on a color image, 
%then saves the output image as a jpg.
%inputs:
%   I:      n x m array of the grayscale values for an image
%   sigma:  Parameters in the Gaussian equation, (1/(2pi*sigma^2))*exp(-(x^2+y^2)/(2*sigma^2)); 
%   tend:   end time
%   dt:     time step
%   str:    string with the name of the image, to be used in the title and
%           saved file.
%
%outputs:
%   I:      The inputted image, after being operated on by heat equation



[m,n] = size(I);

%perform all operations in double precision
I = double(I);
% subplot(121);
% imagesc(I); 
sigma = sqrt(2*tend);


str = sprintf('%s at t=%f', strIn, 0 );
title(str);

A = p;
B = p;
[X,Y] = meshgrid(1:A, 1:B);
G = (1/(2*pi*sigma^2))*exp( -(( X-round(A/2)).^2 + (Y-round(B/2)).^2)/(2*pi*sigma^2));
G = G / sum(sum(G));

I = conv2(I,G,'same');


% %     subplot(122);
% %     imagesc(I);
% %     colormap(gray);
% %     drawnow;
% end
   
    %display the miage
    imagesc(I); 
    colormap(gray);
    str = sprintf('%s at t=%f with sigma=%f and p=%f', strIn, tend,sigma, p );
%     str = sprintf('%s at t=%f with sigma=%f and p=Full Size', strIn, tend,sigma);

    title(str);
    h = figure(1);
    str2 = sprintf('%s%g and p=%g', strIn, tend,p );
%     str2 = sprintf('%s%g and p=Full', strIn, tend );

    
    %save the image 
    saveas(h,str2,'jpg');
    
end

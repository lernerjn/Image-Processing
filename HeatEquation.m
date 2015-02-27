function [ I ] = HeatEquation( I, tend, dt, strIn )
%Author: Jeremy Lerner, Stony Brook University
%This program runs heat equation on an image, then saves the output image
%as a jpg.
%inputs:
%   I:      n x m array of the grayscale values for an image
%   tend:   end time
%   dt:     time step
%   str:    string with the name of the image, to be used in the title and
%           saved file.
%
%outputs:
%   I:      The inputted image, after being operated on by heat equation


dx = 1;
dy = 1;

[m,n] = size(I);

%perform all operations in double precision
I = double(I);
% subplot(121);
imagesc(I); 
colormap(gray);
% Operate on the image with heat equation, using finite differences
for time=dt:dt:tend
    Iold = I;
    
    
    
    % main loop, using a centered finite difference to approximate the
    % second derivatives in heat equation, % dU/dt = d^2U/dx^2 + d^2U/dy^2.
    % Note, this loop is here for readability, the real work is done in
    % Matlab notation below
%     for i=2:m-1
%         for j=2:n-1
%             I(i,j) = Iold(i,j) + dt*((Iold(i+1,j) - 2*Iold(i,j) + Iold(i-1,j)) / dx^2 + (Iold(i,j+1) - 2*Iold(i,j) + Iold(i,j-1))/dy^2); 
%         end
%     end


    %The most efficient way of implementing the for loops, using double
    %colon notation
    I(2:m-1,2:n-1) = Iold(2:m-1,2:n-1) + dt*((Iold(3:m,2:n-1) - 2*Iold(2:m-1,2:n-1) + Iold(1:m-2,2:n-1)) ...
        / dx^2 + (Iold(2:m-1,3:n) - 2*Iold(2:m-1,2:n-1) + Iold(2:m-1,1:n-2))/dy^2);
%     subplot(122);
%     imagesc(I);
%     colormap(gray);
%     drawnow;
end
   
    %display the miage
    imagesc(I); 
    
    colormap(gray);
    str = sprintf('%s at t=%f with dt=%f', strIn, tend,dt );
    title(str);
    h = figure(1);
    str2 = sprintf('%s%g', strIn, tend );
    
    %save the image 
    saveas(h,str2,'jpg');
    
end



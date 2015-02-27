function [ I ] = HeatEquationNonLin( I, tend, dt, strIn )
%Author: Jeremy Lerner, Stony Brook University
%This program runs the non linear heat equation on an image, 
%then saves the output image as a jpg.
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
    
    %homogeneous Neumann boundary conditions, second order accurate
    % that is, the rate of flow at the edges is zero 
%     I(1,:) = Iold(1,:) + 2*(dt/dx^2)*(Iold(2,:) - Iold(1,:));
%     I(end,:) = Iold(end,:) + 2*(dt/dx^2)*(Iold(end,:) - Iold(end-1,:));
%     
%     I(:,1) = Iold(:,1) + 2*(dt/dx^2)*(Iold(:,2) - Iold(:,1));
%     I(:,end) = Iold(:,end) + 2*(dt/dx^2)*(Iold(:,end) - Iold(:,end-1));
    
    
    % main loop, using a centered finite difference to approximate the
    % second derivatives in heat equation, % dU/dt = d^2U/dx^2 + d^2U/dy^2.
    % Note, this loop is here for readability, the real work is done in
    % Matlab notation below
%     for i=2:m-1
%         for j=2:n-1
%             I(i,j) = Iold(i,j) + dt*((Iold(i+1,j) - 2*Iold(i,j) + Iold(i-1,j))*(Iold(i,j+1)-Iold(i,j-1)).^2 / (2*dy*dx^2) ...
%                 + 2*( (Iold(i+1,j+1) - Iold(i-1,j+1) + Iold(i-1,j-1) - Iold(i+1,j-1))*(Iold(i+1,j) - Iold(i-1,j))*(Iold(i,j+1)-Iold(i,j-1)))/(16*(dx^2)*(dy^2)) ...
%                 + (Iold(i,j+1) - 2*Iold(i,j) + Iold(i,j-1))*(Iold(i+1,j)-Iold(i-1,j).^2)/(2*dx*dy^2))/( ((Iold(i+1,j) - Iold(i-1,j))/(2*dx)).^2 + ((Iold(i,j+1) - Iold(i,j-1))/(2*dy)).^2);
%         end
%     end

%     for i=2:m-1
%         for j=2:n-1
%             Uxx = ( Iold(i-1,j) - 2*Iold(i,j) + Iold(i+1,j))/(dx^2);
%             Uyy = ( Iold(i,j-1) - 2*Iold(i,j) + Iold(i,j+1))/(dy^2);
%             Uxy = ( Iold(i+1,j+1) - Iold(i-1,j+1) + Iold(i-1,j-1) - Iold(i+1,j-1)) / (4*dx*dy);
%             Ux  = ( Iold(i+1,j) - Iold(i-1,j)) / (2*dx); 
%             Uy  = ( Iold(i,j+1) - Iold(i,j-1)) / (2*dy); 
%             
%             I(i,j) = Iold(i,j) + dt*( Uxx * Uy^2 - 2*Uxy*Ux*Uy + Uyy*Ux^2)/(Ux^2+Uy^2+1e-14);
%             
%         end
%     end

    %The most efficient way of implementing the for loops, using double
    %colon notation
    
    Uxx = ( Iold(3:m,2:n-1) - 2*Iold(2:m-1,2:n-1) + Iold(1:m-2,2:n-1))./(dx^2);
    Uyy = ( Iold(2:m-1,3:n) - 2*Iold(2:m-1,2:n-1) + Iold(2:m-1,1:n-2))./(dy^2);
    Uxy = ( Iold(3:m,3:n) - Iold(1:m-2,3:n) + Iold(1:m-2,1:n-2) - Iold(3:m,1:n-2))./(4*dx*dy);
    Uy = ( Iold(2:m-1,3:n) - Iold(2:m-1,1:n-2))./(2*dy);
    Ux = ( Iold(3:m,2:n-1) - Iold(1:m-2,2:n-1))./(2*dx);
    

    I(2:m-1,2:n-1) = Iold(2:m-1,2:n-1) + dt.*(Uxx .* Uy.^2 - 2*Uxy.*Ux.*Uy + Uyy.*Ux.^2)./((Ux.^2 + Uy.^2 + 1e-14));
    
%     subplot(122);
%     imagesc(I);
%     colormap(gray);
%     drawnow;
end
   
%     %display the miage
%     imagesc(I); 
%     
%     colormap(gray);
%     str = sprintf('%s at t=%f with dt=%f', strIn, tend,dt );
%     title(str);
%     h = figure(1);
%     str2 = sprintf('%s%g', strIn, tend );
%     
    %save the image 
    %saveas(h,str2,'jpg');
    
end



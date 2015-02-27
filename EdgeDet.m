function [ I ] = EdgeDet( image , save , saveStr, noise, Iinput)
%Author: Jeremy Lerner, Stony Brook University
% Edge Detector based on geometric active contours, or snake based
% segmentation, run on various potential shapes if no input picture is provided.
% Iterate the equation dI/dt = phi *||grad(I)|| *[ div( grad(I) /
% (sqrt(epsilon+||grad(I)||)) + nu] + dotproduct(grad(I),grad(phi))
% until dI/dt = 0. 
%
%   Simple examples: [ I ] = EdgeDet( k , 0 , '', 0, 0); for k=1..12
%
%
% Inputs:  
%       image: an integer between 1 and 12 inclusive to select which image
%       to find the edges of. There are 12 preprogrammed images in this
%       program. 
%
%       save: 1 for saving as an output image, anything else for not
%       
%       saveStr: the string to attach to saved images
%
%       noise: the amount of random noise added to non black pixels,
%       0 adds none, 1 adds very little, 255 is the theoretical maximum,
%       but it likely entirely destroy the image
%
%       I: optional image input, if an image is inputted, the program will
%       ignore image and just work on the inputted image.
% Outputs:
%       I: The zero level set of this output matrix shows the edges of I


%size of I, if I is not provided, m rows, n columns
close all;
m = 200;
n = 200;
nu = 10;
%epsilon as in the equation above, to avoid division by zero or small
%numbers
epsilon = 1;

if nargin == 2
    if save ~= 1
        saveStr = 'blah';
    end
    if save == 1
        saveStr = 'SavedFile';
    end
    noise = 0;
end

if nargin == 3
    noise = 0;
end

if nargin <= 4
    Iinput = 0;
end
    


%dx=dy=1 in the differential equation
dx = 1;
dy = 1;
%time step
dt = 0.05;

%make sure the image is a double
I = 255*ones(m,n);
I = double(I);

%Fill in some shapes in an image. This program will find the edges of the
%shapes, starting from the outside and moving in. If two shapes overlap,
%the program will find the edge of the new overlapped shape, not each shape
%individually
for i = 1:m
    for j =1:n
        

        switch image
            %Weird One
            case 1  
                %first rectangle, top left
                if i > floor(3*m/16) && i < floor(11*m/16) && j > floor(4*n/16) && j < floor(9*n/16)
                    I(i,j) = 0;
                end

                if i > floor(7*m/16) && i < floor(12*m/16) && j > floor(10*i/16) && j < floor(15*i/16)
                     I(i,j) = 0;
                end

                %circle, bottom right
                if  (i - floor(5*m/16))^2 + (j - floor(8*n/16))^2  < 30^2
                    I(i,j) = 0;
                end
                nu = 10;
            %Weird Two
            case 2
                        
                % a rectangle to overlap the circle
                if i > floor(7*m/16) && i < floor(12*m/16) && j > floor(7*n/16) && j < floor(12*n/16)
                    I(i,j) = 0;
                end

                %circle, bottom right
                if  (i - floor(6*m/16))^2 + (j - floor(8*n/16))^2  < 30^2
                    I(i,j) = 0;
                end

                %circle, bottom right
                if  (i - floor(11*m/16))^2 + (j - floor(8*n/16))^2  < 30^2
                    I(i,j) = 0;
                end
                nu = 10;
                
            %Weird Three
            case 3
                if i > floor(7*m/16) && i < floor(12*m/16) && j > floor(8*i/16) && j < floor(15*i/16)
                     I(i,j) = 0;
                end

                if  (i - floor(6*m/16))^2 + (j - floor(8*n/16))^2  < 30^2
                    I(i,j) = 0;
                end

                if i > floor(5*j/16) && i < floor(10*j/16) && j > floor(7*n/16) && j < floor(12*n/16)
                     I(i,j) = 0;
                end
                nu = 10;
                
            %Weird Four
            case 4
                if i > floor(7*m/16) && i < floor(12*m/16) && j > floor(5*i/16) && j < floor(7*i/16)
                     I(i,j) = 0;
                end

                if i > floor(7*m/16) && i < floor(8*m/16) && j > floor(8*i/16) && j < floor(9*i/16)
                     I(i,j) = 0;
                end

                if  (i - floor(6*m/16))^2 + (j - floor(8*n/16))^2  < 30^2
                    I(i,j) = 0;
                end

                if i > floor(5*j/16) && i < floor(10*j/16) && j > floor(7*n/16) && j < floor(12*n/16)
                     I(i,j) = 0;
                end

                if i > floor(6.5*m/16) && i < floor(7.5*m/16) && j > floor(3*n/16) && j < floor(7*n/16)
                    I(i,j) = 0;
                end

                if  (i - floor(12*m/16))^2 + (j - floor(4*n/16))^2  < 10^2
                    I(i,j) = 0;
                end
                nu = 10;
            
            %single circle, bottom right
            case 5
                if (i - floor(9*m/16))^2 + (j - floor(9*n/16))^2  < 45^2
                    I(i,j) = 0;
                end
                nu = 10;
            %single rectangle, top left
            case 6
                if i > floor(4*m/16) && i < floor(11*m/16) && j > floor(4*n/16) && j < floor(10*n/16)
                    I(i,j) = 0;
                end
                nu = 10;
                
                
            %grayscale circle that fades
            case 7
                if (i - floor(10*m/20))^2 + (j - floor(9*n/20))^2  < 70^2 
                    dist = (255/70)*sqrt( (i - 10*m/20)^2 + (j - 9*n/20)^2);
                    I(i,j) = dist;
                end
                
                if (i - floor(10*m/20))^2 + (j - floor(9*n/20))^2  < 20^2
                    I(i,j) = 0;
                end
                
                nu = 10;
                
            %gradient shade
            case 8
                if i > 3*m/20 && i < 17*m/20 && j > 5*n/20 && j <= 7*n/20
                    I(i,j) = 0;
                end
                if i > 3*m/20 && i < 17*m/20 && j > 7*n/20 && j < 11*n/20
                    I(i,j) = j*255/(11*n/20);
                end
                nu = 10;
                 

            %thin shell, 2 pixels wide
            case 9
                if i > 6.8*m/20 && i < 11.1*m/20 && j > 4.8*n/20 && j < 5.1*n/20
                    I(i,j) = 0;
                end
                if i > 7*m/20 && i < 11.1*m/20 && j > 10.8*n/20 && j < 11.1*n/20
                    I(i,j) = 0;
                end
                if i > 6.8*m/20 && i < 7.1*m/20 && j > 5*n/20 && j < 11.1*n/20
                    I(i,j) = 0;
                end
                
                 if i > 10.8*m/20 && i < 11.1*m/20 && j >= 5*n/20 && j < 11.1*n/20
                    I(i,j) = 0;
                 end
                 nu = 10;
                 
            %thicker shell, 5 pixels wide
            case 10
                if i > 6.6*m/20 && i < 10.1*m/20 && j > 4.6*n/20 && j < 5.1*n/20
                    I(i,j) = 0;
                end
                if i > 7*m/20 && i < 10.1*m/20 && j > 9.6*n/20 && j < 10.1*n/20
                    I(i,j) = 0;
                end
                if i > 6.6*m/20 && i < 7.1*m/20 && j > 5*n/20 && j < 10.1*n/20
                    I(i,j) = 0;
                end
                
                 if i > 9.6*m/20 && i < 10.1*m/20 && j >= 5*n/20 && j < 10.1*n/20
                    I(i,j) = 0;
                 end
                 nu = 10;
%             %very thin shell
%             case 10
%                 if i > 6.95*m/20 && i < 10.1*m/20 && j > 4.95*n/20 && j < 5.1*n/20
%                     I(i,j) = 0;
%                 end
%                 if i > 7*m/20 && i < 10.1*m/20 && j > 9.95*n/20 && j < 10.1*n/20
%                     I(i,j) = 0;
%                 end
%                 if i > 6.95*m/20 && i < 7.1*m/20 && j > 5*n/20 && j < 10.1*n/20
%                     I(i,j) = 0;
%                 end
%                 
%                  if i > 9.95*m/20 && i < 10.1*m/20 && j >= 5*n/20 && j < 10.1*n/20
%                     I(i,j) = 0;
%                  end
%                 nu = 10;

                case 11  
                    %first rectangle, top left
                    if i > floor(3*m/16) && i < floor(11*m/16) && j > floor(3*n/16) && j < floor(9*n/16)
                        I(i,j) = 0;
                    end

                    %circle, bottom right
                    if  (i - floor(9*m/16))^2 + (j - floor(12*n/16))^2  < 25^2
                        I(i,j) = 0;
                    end
                    nu = 10;
                
                case 12  
                    %first rectangle, top left
                    if  (i - floor(5*m/16))^2 + 5*(j - floor(5*n/16))^2  < 25^2
                        I(i,j) = 0;
                    end

                    %circle, bottom right
                    if  (i - floor(9*m/16))^2 + (j - floor(12*n/16))^2  < 25^2
                        I(i,j) = 0;
                    end
                    nu = 10;
        end
               
    end
end
%use the user inputted image
if Iinput ~= 0
    I = Iinput;
    
end



%Add random noise to the image if the user requests it, except do not add
%noise to black pixels, because that causes numerical problems
if ( noise ~= 0)
    n2 = noise;

    for i = 1:m
        for j = 1:n
            %make sure to not add noise that makes the pixels out of range
            if I(i,j) < 255 - n2 
                I(i,j) = I(i,j) + n2*rand(1);

            elseif I(i,j) > n2
                I(i,j) = I(i,j) - n2*rand(1);
            end
            
        end
    end

end

% Optional (removed for now): smoothing with convolution, 
% the size of the convolution matrix makes a
% large difference in the smoothing, can be uncommented
%    G = fspecial('gaussian',[5,5],1);
% % %   # Filter it
%   I = imfilter(I,G,'same');

% Smoothing with nonlinear heat equation
% [ I ] = double(abs(HeatEquationNonLin( I, 2, 0.1 )));

%save the original image
Iorig = I;



%plot the original image
subplot(131);
imagesc(I); 
title('The Original Image');
colormap(gray);
phi = double(ones(m,n));
% 

% Note: There are many methods of calculating phi numerically, for this
% project, I used second order centered differences in the end.
% But many options are programmed in, just commented out. Both higher
% and lower order accuracy are available
% 
% % phi = 1/(1+ || grad(I)||^2) , the "edginess" of each pixel, using a
% CENTERED difference. Because the gradient is large near edges, the
% function phi is small near edges and nearly one away from edges.
for i = 2:size(I,1)-1
    for j = 2:size(I,2)-1
        phi(i,j) = 1 / ( 1 + ((I(i+1,j) - I(i-1,j))./(2*dx)).^2 + ((I(i,j+1) - I(i,j-1))./(2*dy)).^2);
    end
end
% % 
% % % phi = 1/(1+ || grad(I)||^2) , the "edginess" of each pixel, using a
% % Fourth order CENTERED difference. Because the gradient is large near edges, the
% % function phi is small near edges and nearly one away from edges.
% for i = 3:size(I,1)-2
%     for j = 3:size(I,2)-2
%         phi(i,j) = 1 / ( 1 + ((-I(i+2,j) + 8*I(i+1,j) - 8*I(i-1,j) + I(i-2,j))./(12*dx)).^2 ... 
%             + ((-I(i,j+2) + 8*I(i,j+1) - 8*I(i,j-1) + I(i,j-2))./(12*dy)).^2);
%     end
% end
% % 
% % % phi = 1/(1+ || grad(I)||^2) , the "edginess" of each pixel, using a
% % FORWARD difference. Because the gradient is large near edges, the
% % function phi is small near edges and nearly one away from edges.
% for i = 2:size(I,1)-1
%     for j = 2:size(I,2)-1
%         phi(i,j) = 1 / ( 1 + ((I(i+1,j) - I(i,j))./(dx)).^2 + ((I(i,j+1) - I(i,j))./(dy)).^2);
%     end
% end

% % % phi = 1/(1+ || grad(I)||^2) , the "edginess" of each pixel, using a
% % BACKWARD difference. Because the gradient is large near edges, the
% % function phi is small near edges and nearly one away from edges.
% for i = 2:size(I,1)-1
%     for j = 2:size(I,2)-1
%         phi(i,j) = 1 / ( 1 + ((I(i,j) - I(i-1,j))./(dx)).^2 + ((I(i,j) - I(i,j-1))./(dy)).^2);
%     end
% end


% % % phi = 1/(1+ || grad(I)||^2) , the "edginess" of each pixel, using a
% % second order FORWARD difference. Because the gradient is large near edges, the
% % function phi is small near edges and nearly one away from edges.
% for i = 2:size(I,1)-2
%     for j = 2:size(I,2)-2
%         phi(i,j) = 1 / ( 1 + ((-3*I(i,j) + 4*I(i+1,j) - I(i+2,j))./(2*dx)).^2 ...
%             + ((-3*I(i,j) +4 *I(i,j+1)-I(i,j+2))./(2*dy)).^2);
%     end
% end


% % % phi = 1/(1+ || grad(I)||^2) , the "edginess" of each pixel, using a
% % NINTH order CENTERED difference. Because the gradient is large near edges, the
% % function phi is small near edges and nearly one away from edges.
% for i = 5:size(I,1)-4
%     for j = 5:size(I,2)-4
%         phi(i,j) = 1 /( 1 + ((3*I(i-4,j) - 32*I(i-3,j) + 168*I(i-2,j) - 672*I(i-1,j) ...
%             + 672*I(i+1,j) - 168*I(i+2,j) + 32*I(i+3,j) - 3*I(i+4,j) )./(840*dx)).^2 ...
%             + ((3*I(i,j-4) - 32*I(i,j-3) + 168*I(i,j-2) - 672*I(i,j-1) ...
%             + 672*I(i,j+1) - 168*I(i,j+2) + 32*I(i,j+3) - 3*I(i,j+4) )./(840*dx)).^2);
%     end
% end

% % % phi = 1/(1+ || grad(I)||^2) , the "edginess" of each pixel, using a
% % SEVENTH order CENTERED difference. Because the gradient is large near edges, the
% % function phi is small near edges and nearly one away from edges.
% for i = 4:size(I,1)-3
%     for j = 4:size(I,2)-3
%         phi(i,j) = 1 /( 1 + ((-I(i-3,j) + 9*I(i-2,j) -45*I(i-1,j) + 45*I(i+1,j) - 9*I(i+2,j) + I(i+3,j) )./(60*dx)).^2 ...
%             +((-I(i,j-3) + 9*I(i,j-2) -45*I(i,j-1) + 45*I(i,j+1) - 9*I(i,j+2) + I(i,j+3) )./(60*dx)).^2);
%     end
% end



%The snake originally, a circle with a center at the middle of the image.
[X,Y] = meshgrid(1:size(I,1), 1:size(I,2));
InitCond = double((X-m/2).^2 + (Y-n/2).^2 - 8000)/(1000) ;
%a different possible initial condition
% InitCond = abs(X-m/2)+abs(Y-n/2)-100;
% InitCond = double(InitCond)/10000;

% dphi/dx, the derivative in the horizontal direction
phix= diffxy(1,phi,2,1);
% dphi/dy, the derivative in the vertical direction
phiy= diffxy(1,phi,1,1);

% Used for upwinding the third term in the PDE
aplusx = max(-phix,0);
aminusx = min(-phix,0);
aplusy = max(-phiy,0);
aminusy = min(-phiy,0);


%Plot the original snake
I = InitCond;
subplot(132);
contour(I, [0,0]); 
title('The original snake');
set(gca,'YDir','Reverse');
colormap(gray);


Iold = 255*ones(m,n);

iterations = 0;
Iolder =  ones(10,10);
Ioldest = ones(10,10);
        

%if I has stopped changing, dI/dt = 0, stop iterating
while norm(I - Iold) > 1e-3
%   rebuilding, every five iterations
    if mod(iterations,5) == 0
        Ioldest = Iold;
        
        I = double(redistanceV(I));
        
    end
    
    Iold = double(I);
   
    %Use a second order accurate estimate of the first and second
    %derivatives, see comments below (second order accurate in the interior
    %and at the boundaries)
    Ux = diffxy(1,Iold,2,1);
    Uy = diffxy(1,Iold,1,1);
    Uxx = diffxy(1,Iold,2,2);
    Uyy = diffxy(1,Iold,1,2);
    Uxy = diffxy(1,Ux,1,1);
    
    uxminus = zeros(m,n);
    uxplus = zeros(m,n);
   
    % backward difference, used in upwinding in the PDE
    uxminus(1:m,2:n) = Iold(1:m,2:n) - Iold(1:m, 1:n-1);
    uxminus(:,1) = ((-3*Iold(:,1) + 4*Iold(:,2) - Iold(:,3) )/ 2)  ;
    % forward difference
    uxplus(1:m,1:n-1) = Iold(1:m,2:n) - Iold(1:m, 1:n-1);
    uxplus(:,n) = ((-3*Iold(:,n-2) + 4*Iold(:,n-1) - Iold(:,n) )/ 2)  ;
    
    
    uyminus = zeros(m,n);
    uyplus = zeros(m,n);

    % backward difference, used in upwinding in the PDE
    uyminus(2:m,1:n) = Iold(2:m,1:n) - Iold(1:m-1, 1:n);
    uyminus(1,:) = ((-3*Iold(1,:) + 4*Iold(2,:) - Iold(3,:) )/ 2)  ;
    % forward difference
    uyplus(1:m-1,1:n) = Iold(2:m,1:n) - Iold(1:m-1, 1:n);
    uyplus(m,:) = ((-3*Iold(:,n-2) + 4*Iold(:,n-1) - Iold(:,n) )/ 2)  ;
    
    % The main part of the loop, using finite differences for 
    %   dI/dt = phi *||grad(I)|| *[ div( grad(I) /
    %   (sqrt(epsilon+||grad(I)||))] phi*||grad(I)||*nu + dotproduct(grad(I),grad(phi))
    % with second order centered differences for the first term, and
    % entropy conservating upwinding solutions for the second and third
    % terms
    I = Iold + dt*phi.*((Ux.^2 + Uy.^2).^(1/2)).*(((1./(epsilon+Ux.^2+Uy.^2).^(3/2)).*...
   (epsilon.*(Uxx + Uyy) + (Ux.^2).*Uyy + (Uy.^2).*Uxx - 2.*(Ux).*(Uy).*Uxy))) ...
   + dt*nu.*phi.*sqrt((max(uxplus,0)).^2 + min(uxminus,0).^2 + (max(uyplus,0)).^2 + (min(uyminus,0)).^2) ...
   - dt.*(aplusx.*uxminus + aminusx.*uxplus + aplusy.*uyminus + aminusy.*uyplus ) ;



    % Homogeneous Neumann boundary conditions, dI/dt = 0 at the edge
    I(1,:) = Iold(2,:);
    I(m,:) = Iold(m-1,:);
    I(:,1) = Iold(:,2);
    I(:,n) = Iold(:,n-1);


    %plot the current progress
    subplot(133);
    contour(I, [0,0]);
    set(gca,'YDir','Reverse');
    colormap(gray); 
    title(['The snake after ' num2str(iterations) ' iterations']);
    drawnow
%     
% % %     The current state of the contour, useful for debugging
%     figure(2);
%     surf(I);
%     title(['The snake after ' num2str(iterations) ' iterations']);
%     drawnow
% %     
    iterations = iterations + 1;
    tol = 1e-5;
%     %if the contour has stopped changing, relative to the last iteration, 
%     % the one before that or the one before redistancing
    if ( size(contourc(I, [0,0])) == size(contourc(Iold, [0,0]))) 
        norm(norm(contourc(I, [0,0]) - contourc(Iold, [0,0])));
        if( norm(norm(contourc(I, [0,0]) - contourc(Iold, [0,0]))) < tol)
            break
        end
    end
    if ( size(contourc(I, [0,0])) == size(contourc(Iolder, [0,0]))) 
        norm(norm(contourc(I, [0,0]) - contourc(Iolder, [0,0])));
        if (norm(norm(contourc(I, [0,0]) - contourc(Iolder, [0,0]))) < tol)
            break
        end
    end 
    
    if ( size(contourc(Iold, [0,0])) == size(contourc(Iolder, [0,0]))) 
        norm(norm(contourc(Iold, [0,0]) - contourc(Iolder, [0,0])));
        if (norm(norm(contourc(Iold, [0,0]) - contourc(Iolder, [0,0]))) < tol)
            break
        end
    end
    
    if ( size(contourc(Ioldest, [0,0])) == size(contourc(I, [0,0]))) 
        norm(norm(contourc(Ioldest, [0,0]) - contourc(I, [0,0])));
        if (norm(norm(contourc(Ioldest, [0,0]) - contourc(I, [0,0]))) < tol)
            break
        end
    end
   
    
    Iolder = Iold;
end


if ( save == 1)

    h = figure(1);
    str = sprintf('SqAndCirc%g' , iterations );

    %save the image 
    saveas(h,str,'jpg');
end    

figure(2)
colormap(gray);
imagesc(Iorig);

hold on;
[~,h] = contour(I, [0,0], 'g');
set(h,'linewidth',2);
str2 = sprintf('The Final Contour and Original Image, after %g iterations', iterations);
title(str2);

if ( 1 && save == 1 )
    h = figure(2);
%     str3 = sprintf('ThinShell_Final%g' , nu );
    saveas(h,saveStr,'jpg');
    
end

if (0 &&  save == 1)
    figure(3)
    surf(I);
    hold on;
    P = contourc(I,[0,0]);
    line(P(1,2:end), P(2,2:end),zeros(size(P,2)-1),'linewidth',2,'color','r')
    h = figure(3);
    str3 = sprintf('CircSurf%g' , iterations );
    saveas(h,str3,'jpg');
end



end

%Deceptively named, this function actually rebuilds PHI
% Vadim Ratner's redistancing/rebuilding function
function u = redistanceV(Phi)
    mask_in = (Phi<=0); %inside the contour
    mask_out = 1-mask_in; %outside the contour
    u = (double(bwdist(mask_in)) +(1-double(bwdist(mask_out))).*mask_in);
    u = double(u);
end



%Smoothing function, if necessary
function [ I ] = HeatEquationNonLin( I, tend, dt )
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
    
    Uxx = double(( Iold(3:m,2:n-1) - 2*Iold(2:m-1,2:n-1) + Iold(1:m-2,2:n-1))./(dx^2));
    Uyy = double(( Iold(2:m-1,3:n) - 2*Iold(2:m-1,2:n-1) + Iold(2:m-1,1:n-2))./(dy^2));
    Uxy = double(( Iold(3:m,3:n) - Iold(1:m-2,3:n) + Iold(1:m-2,1:n-2) - Iold(3:m,1:n-2))./(4*dx*dy));
    Uy = double(( Iold(2:m-1,3:n) - Iold(2:m-1,1:n-2))./(2*dy));
    Ux = double(( Iold(3:m,2:n-1) - Iold(1:m-2,2:n-1))./(2*dx));
    

    I(2:m-1,2:n-1) = Iold(2:m-1,2:n-1) + double(dt.*(Uxx .* Uy.^2 - 2*Uxy.*Ux.*Uy + Uyy.*Ux.^2)./((Ux.^2 + Uy.^2 + 1e-14)));

end
   
   
end



%Imported and borrowed numerical derivative function
function dy = diffxy(x,y,varargin)
% DIFFXY - accurate numerical derivative/differentiation of Y w.r.t X. 
%
%   DY = DIFFXY(X,Y) returns the derivative of Y with respect to X using a 
%        pseudo second-order accurate method. DY has the same size as Y.
%   DY = DIFFXY(X,Y,DIM) returns the derivative along the DIM-th dimension
%        of Y. The default is differentiation along the first 
%        non-singleton dimension of Y.
%   DY = DIFFXY(X,Y,DIM,N) returns the N-th derivative of Y w.r.t. X.
%        The default is 1.
%
%   Y may be an array of any dimension.
%   X can be any of the following:
%       - array X with size(X) equal to size(Y)
%       - vector X with length(X) equal to size(Y,DIM)
%       - scalar X denotes the spacing increment
%   DIM and N are both integers, with 1<=DIM<=ndims(Y)
%
%   DIFFXY has been developed especially to handle unequally spaced data,
%   and features accurate treatment for end-points.
%
%   Example: 
%   % Data with equal spacing
%     x = linspace(-1,2,20);
%     y = exp(x);
% 
%     dy = diffxy(x,y);
%     dy2 = diffxy(x,dy);  % Or, could use >> dy2 = diffxy(x,y,[],2);
%     figure('Color','white')
%     plot(x,(y-dy)./y,'b*',x,(y-dy2)./y,'b^')
%
%     Dy = gradient(y)./gradient(x);
%     Dy2 = gradient(Dy)./gradient(x);
%     hold on
%     plot(x,(y-Dy)./y,'r*',x,(y-Dy2)./y,'r^')
%     title('Relative error in derivative approximation')
%     legend('diffxy: dy/dx','diffxy: d^2y/dx^2',...
%            'gradient: dy/dx','gradient: d^2y/dx^2')
%
%   Example: 
%   % Data with unequal spacing. 
%     x = 3*sort(rand(20,1))-1;
%     % Run the example above from y = exp(x)
%
%   See also DIFF, GRADIENT
%        and DERIVATIVE on the File Exchange

% for Matlab (should work for most versions)
% version 1.0 (Nov 2010)
% (c) Darren Rowland
% email: darrenjrowland@hotmail.com
%
% Keywords: derivative, differentiation

[h,dy,N,perm] = parse_inputs(x,y,varargin);
if isempty(dy)
    return
end
n = size(h,1);
i1 = 1:n-1;
i2 = 2:n;

for iter = 1:N
    v = diff(dy)./h;
    if n>1
        dy(i2,:) = (h(i1,:).*v(i2,:)+h(i2,:).*v(i1,:))./(h(i1,:)+h(i2,:));
        dy(1,:) = 2*v(1,:) - dy(2,:);
        dy(n+1,:) = 2*v(n,:) - dy(n,:);
    else
        dy(1,:) = v(1,:);
        dy(n+1,:) = dy(1,:);
    end
end

% Un-permute the derivative array to match y
dy = ipermute(dy,perm);

end
%%% Begin local functions %%%
function [h,dy,N,perm] = parse_inputs(x,y,v)

numvarargs = length(v);
if numvarargs > 2
    error('diffxy:TooManyInputs', ...
        'requires at most 2 optional inputs');
end

h = [];
N = [];
perm = [];

% derivative along first non-singleton dimension by default
dim = find(size(y)>1);
% Return if dim is empty
if isempty(dim)
    dy = [];
    return
end
dim = dim(1);

% Set defaults for optional arguments
optargs = {dim 1};
newVals = ~cellfun('isempty', v);
optargs(newVals) = v(newVals);
[dim, N] = optargs{:};

% Error check on inputs
if dim<1 || dim>ndims(y) || dim~=fix(dim) || ~isreal(dim)
    error('diffxy:InvalidOptionalArg',...
        'dim must be specified as a non-negative integer')
end
if N~=fix(N) || ~isreal(N)
    error('diffxy:InvalidOptionalArg',...
        'N must be an integer')
end

% permutation which will bring the target dimension to the front
perm = 1:length(size(y));
perm(dim) = [];
perm = [dim perm];
dy = permute(y,perm);


if length(x)==1  % Scalar expansion to match size of diff(dy,[],1)
    sizeh = size(dy);
    sizeh(1) = sizeh(1) - 1;
    h = repmat(x,sizeh);
elseif ndims(x)==2 && any(size(x)==1) % Vector x expansion
    if length(x)~=size(dy,1)
        error('diffxy:MismatchedXandY',...
            'length of vector x must match size(y,dim)')
    end
    x = x(:);
    sizeh = size(dy);
    sizeh(1) = 1;
    h = repmat(diff(x),sizeh);
else
    if size(y) ~= size(x)
        error('diffxy:MismatchedXandY',...
            'mismatched sizes of arrays x and y');
    end
    % Permute x as for y, then diff
    h = diff(permute(x,perm),[],1);
end
end



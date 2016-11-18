function [I] = Create_Seg_Image(image, noise, m, n)

%make sure the image is a double
I = 255*ones(m,n);
I = double(I);
%Fill in some shapes in an image. This program will find the edges of the
%shapes, starting from the outside and moving in. If two shapes overlap,
%the program will find the edge of the new overlapped shape, not each shape
%individually

% Optional: fill in the bottom and top half of the image separately
% I( 1 : end/2 - 1, :) = 255;% I( 1 : end/2 - 1, :) - 50;
% I( end/2 : end, :) = 0; %I( end/2 : end, :) + 20;

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
            
            %single circle, bottom right
            case 5
                if (i - floor(8*m/16))^2 + (j - floor(8*n/16))^2  < 75^2 %45: old value
                    I(i,j) = 127;
                end
            %single rectangle, top left
            case 6
                if i > floor(4*m/16) && i < floor(11*m/16) && j > floor(4*n/16) && j < floor(10*n/16)
                    I(i,j) = 0;
                end

                
                
            %grayscale circle that fades
            case 7
                if (i - floor(10*m/20))^2 + (j - floor(9*n/20))^2  < 70^2 
                    dist = (255/70)*sqrt( (i - 10*m/20)^2 + (j - 9*n/20)^2);
                    I(i,j) = dist;
                end
                
                if (i - floor(10*m/20))^2 + (j - floor(9*n/20))^2  < 20^2
                    I(i,j) = 0;
                end
                

            %gradient shade
            case 8
                if i > 3*m/20 && i < 17*m/20 && j > 5*n/20 && j <= 7*n/20
                    I(i,j) = 0;
                end
                if i > 3*m/20 && i < 17*m/20 && j > 7*n/20 && j < 11*n/20
                    I(i,j) = j*255/(11*n/20);
                end

                 

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
                
                case 12  
                    %first circle, top left
                    if  (i - floor(5*m/16))^2 + 5*(j - floor(5*n/16))^2  < 45^2
                        I(i,j) = 0;
                    end

                    %circle, bottom right
                    if  (i - floor(9*m/16))^2 + (j - floor(12*n/16))^2  < 45^2
                        I(i,j) = 0;
                    end
        end
               
    end
end



% I(I==0) = 100;
pm = rand(size(I));
pm(pm<0.5) = -1;
pm(pm>=0.5) = 1;

pm_2 = rand(size(I));
pm_2(pm_2<0.5) = -1;
pm_2(pm_2>=0.5) = 1;

if ( noise ~= 0)

    for i = 1:m
        for j = 1:n
            
            if I(i,j) < (255 - noise)  && I(i,j) > noise
                I(i,j) = I(i,j) + pm(i,j)*rand(1) +pm_2(i,j)*(noise/2)*rand(1);

            elseif I(i,j) <= noise
                I(i,j) = I(i,j) + noise*rand(1);
            else
                I(i,j) = I(i,j) - noise*rand(1);
            end
            
        end
    end

end




%Add random noise to the image if the user requests it, except do not add
%noise to black pixels, because that causes numerical problems
% if ( noise ~= 0)
%     n2 = noise;
% 
%     for i = 1:m
%         for j = 1:n
%             
%             if I(i,j) < 255 - n2  %&& I(i,j) > 0
%                 I(i,j) = I(i,j) + n2*rand(1);
% 
%             elseif I(i,j) > n2
%                 I(i,j) = I(i,j) - n2*rand(1);
%             end
%             
%         end
%     end
% 
% end
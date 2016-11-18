function phi = CV_segmentation(Img)

% CV_bw: Chan Vese segmentation
%                   for black and white (i.e. single channel) images
% inputs:
%               Img: the image to be segmented
% outputs:
%               phi: the 3d level set function where the zero level set
%               represent the 2d image segmentation
% CreateImage4;
% Img_2D   = double(image_created);
% Img_2D  = (double(rgb2gray( Img_2D )   ));

% Img_path = 'Images/Shapes/s_all_water1.png';
% Img_path = 'Images/Shapes/c2_1.png';
% Img  = (double(rgb2gray(imread(Img_path))));


image_code = 3;
noise = 145;

n1 = size(Img,1);
n2 = size(Img,2) ;
n = n1;
m = n2;

% Img = Create_Seg_Image(image_code, noise, m, n);
            

[rr, cc] = meshgrid(1:size(Img,2), 1:size(Img,1));
C = sqrt((rr-size(Img,2)/2).^2+(cc-size(Img,1)/2).^2)<= 0.75*max(size(Img))/2;

% C = sqrt((rr-size(Img_2D,2)/2).^2+(cc-size(Img_2D,1)/2).^2)<= 50;

level_set_offset = 2;

phi = bwdist(C)-(bwdist(1-C)) -0.5;
% mean_in = mean(Img(phi<=-level_set_offset));
% mean_out = mean(Img(phi >level_set_offset));

mean_out = sum(sum(Img.*Heaviside(phi)))/(sum(sum(Heaviside(phi)))+1e-14);
mean_in = sum(sum(Img.*(1-Heaviside(phi))))/(sum(sum((1-Heaviside(phi))))+1e-14);
    
  
%% Vectorized version of the Bhattacharyya Segmentation

num_iters = 1000;
phi_old = zeros(size(phi));
iterations = 1;
dt = 0.1;
alpha = 1;

[m,n] = size(phi);
dx = 1;
dy = 1;
Ux = zeros(size(phi));
Uy = zeros(size(phi));
Uxy = zeros(size(phi));
Uxx = zeros(size(phi));
Uyy = zeros(size(phi));

while norm(phi-phi_old) > 1e-1 && iterations < num_iters
    
%     if  mod(iterations,1) == 0 || iterations == 1
%         offset_phi = phi + level_set_offset;
%         figure(1)
%         imshow(uint8(Img),'InitialMagnification',65); hold on;
%         [c,h] = contour(offset_phi,[1,1],'y'); set(h,'LineWidth',4); hold off;
%         set(gcf,'Color','w'); axis off; drawnow;
%         title(strcat('mean value in:' , num2str(mean_in), ', mean value out: ' , num2str(mean_out)))
%         
% %         figure(2)
% %         surf(phi)
%         
%     end


    % Curvature of phi
    Uxx(2:m-1,2:n-1) = ( phi(3:m,2:n-1) - 2*phi(2:m-1,2:n-1) + phi(1:m-2,2:n-1))./(dx^2);
    Uyy(2:m-1,2:n-1) = ( phi(2:m-1,3:n) - 2*phi(2:m-1,2:n-1) + phi(2:m-1,1:n-2))./(dy^2);
    Uxy(2:m-1,2:n-1) = ( phi(3:m,3:n) - phi(1:m-2,3:n) + phi(1:m-2,1:n-2) - phi(3:m,1:n-2))./(4*dx*dy);
    Uy(2:m-1,2:n-1) = ( phi(2:m-1,3:n) - phi(2:m-1,1:n-2))./(2*dy);
    Ux(2:m-1,2:n-1) = ( phi(3:m,2:n-1) - phi(1:m-2,2:n-1))./(2*dx);
    curvature = (Uxx .* Uy.^2 - 2*Uxy.*Ux.*Uy + Uyy.*Ux.^2)./((Ux.^2 + Uy.^2 + 1e-14).^(2/2));

    
%     mean_in = mean(Img(phi<=-level_set_offset));
%     mean_out = mean(Img(phi >level_set_offset));

    mean_out = sum(sum(Img.*Heaviside(phi)))/(sum(sum(Heaviside(phi)))+1e-14);
    mean_in = sum(sum(Img.*(1-Heaviside(phi))))/(sum(sum((1-Heaviside(phi))))+1e-14);
    
    
     

    phi_old = phi;
    phi = phi - Delta(phi)*dt.*( ( (double(Img) - mean_out).^2 - (double(Img)-mean_in).^2 ).*sqrt(Ux.^2+Uy.^2) - alpha*curvature);%/max(max(abs(curvature)));
     
    mask_in = (phi <= 0);
    mask_out = 1-mask_in;
    phi = double(double(bwdist(mask_in)) + (1-double(bwdist(mask_out)))) - 0.5;
    
    
    iterations = iterations + 1;

end

        offset_phi = phi + level_set_offset;
        figure(1)
        imshow(uint8(Img),'InitialMagnification',65); hold on;
        [c,h] = contour(offset_phi,[1,1],'y'); set(h,'LineWidth',4); hold off;
        set(gcf,'Color','w'); axis off; drawnow;
        title('Final Segmentation');
%         title(strcat('mean value ins: ' , num2str(mean_in), ', mean value out: ' , num2str(mean_out)))
        

%% Allen Meeting
% integral((I-J)^2)+ integral( || grad(J)||^2), J approximates I, J is
% smooth
% Chan Vese:
% integral( (I-mu_out)^2, over omega_out)) + integral( (I-mu_un)^2, over omega_in))
%
% Variance: E(x^2 - (E(x))^2) = variance = E(x^2) - (E(x))^2 = E(x^2) when mean is zero

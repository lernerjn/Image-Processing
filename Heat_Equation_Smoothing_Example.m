load('Tovi_BW');
load('Tovi_RGB');
I = Tovi_BW;
I2 = Tovi_RGB;
[ ISmoothed ] = HeatEquationNonLin( I, 25, 0.2, 'Tovi' );
[ I2Smoothed ] = HeatEquationNonLinRGB( I2, 25, 0.2, 'Tovi' );


subplot(2,2,1)
imagesc(I);
title('The Original Black and White Image');
subplot(2,2,2)
imagesc(ISmoothed);
title('The Smoothed Black and White Image');
colormap(gray);

subplot(2,2,3)
imagesc(I2);
title('The Original Color Image');
subplot(2,2,4)
imagesc(I2Smoothed);
title('The Smoothed Color Image');

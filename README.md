# Image-Processing
For an example of smoothing using heat equation, run Heat_Equation_Smoothing_Example.m, 
which will automatically display the results.

For examples of edge detection using (edge-based) geodesic active contours (which will automatically display the results), run:  

	[I] = EdgeDet( k , 0, ‘’, 0, 0); % for k = 1 through 12
	
For examples of edge detection using region-based segmentation, run:

	% Optional line to create a black and white image
	% Note that k = 1 through 12 where noise adds randomized noise to the image, 0 <= noise < 255 
	% and the image is of size m x n.
	Img = Create_Seg_Image(k, noise, m, n)
	% Run the actual segmentation on Img
	phi = CV_segmentation(Img);

Background on all codes is provided in the Analysis folder.


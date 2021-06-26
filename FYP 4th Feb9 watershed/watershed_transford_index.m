%% Method of Gabor Filter

lambda  = 5.48;
%for patient5.png, lambda is 6
%for patient1.jpg, lambda is 5.48 and se1 is 3, bwarea = 150, thres = 0.65
%for nodule.jpg/patient4 lambda is 5.48, otherwise 6
%for ct/patient1.jpg lambda 5.48 further improves the result and makes the nodule
%smaller. 
% for patient6.jpg, lambda is 5.48 for patient6.jpg, 4 becomes 3 -se1 in
% line 282, in line 182 = threshold is 0.82
theta   = 0;
bw      = 3;
psi     = [0 0];
gamma   = 2;
N       = 4;
%% Ct scan input image
 img_in = imread('F:\Local Disk D\Study\Final Year Project\Final Project\FYP 4th Feb9 watershed\Dataset\patient1.jpg');
 img_in = imresize(img_in,[512 512]);
 if size(img_in,3)==3
    img_in = rgb2gray(img_in);
end 
%img_in = b;
%img_in = double(dicomread('a.dcm'));
%img_in(:,:,2:3) = [];
img_out = zeros(size(img_in,1), size(img_in,2), N);
for n=1:N
    gb = gabor_fn(bw,gamma,psi(1),lambda,theta)...
        + gabor_fn(bw,gamma,psi(2),lambda,theta);
    img_out(:,:,n) = imfilter(img_in, gb, 'symmetric');
    theta = theta + 2*pi/N;
end
figure(1);
imshow(img_in);
title('input image');
figure(2);
img_out_disp = sum(abs(img_out).^2, 3).^0.5;
img_out_disp = img_out_disp./max(img_out_disp(:));
imshow(img_out_disp);
I = img_out_disp;
%a=imcrop(img_out_disp,[65 120 390 300]);
%a=img_out_disp(120:400,65:466);
%a=img_out_disp(150:420,45:470);
title('gabor output, L-2 super-imposed, normalized');
%% Step 1: Read in the Color Image and Convert it to Grayscale

%I=a;
%I = img_out_disp;

%% Step 2: Use the Gradient Magnitude as the Segmentation Function
% Use the Sobel edge masks, |imfilter|, and some simple arithmetic to
% compute the gradient magnitude.  The gradient is high at the borders of
% the objects and low (mostly) inside the objects.

%preprocessing using sobel filter

hy = fspecial('sobel');
hx = hy';
Iy = imfilter(double(I), hy, 'replicate');
Ix = imfilter(double(I), hx, 'replicate');
gradmag = sqrt(Ix.^2 + Iy.^2);
figure
imshow(gradmag,[]), title('Gradient magnitude (gradmag)')

%% 
% Can you segment the image by using the watershed transform directly on
% the gradient magnitude?


L = watershed(gradmag);
Lrgb = label2rgb(L);
figure, imshow(Lrgb), title('Watershed transform of gradient magnitude (Lrgb)')

%%
% No.  Without additional preprocessing such as the marker computations below,
% using the watershed transform directly often results in
% "oversegmentation."


%% Step 3: Mark the Foreground Object
% A variety of procedures could be applied here to find the foreground
% markers, which must be connected blobs of pixels inside each of the
% foreground objects.  In this example you'll use morphological
% techniques called "opening-by-reconstruction" and
% "closing-by-reconstruction" to "clean" up the image.  These
% operations will create flat maxima inside each object that can be
% located using |imregionalmax|.


%%
% Opening is an erosion followed by a dilation, while
% opening-by-reconstruction is an erosion followed by a morphological
% reconstruction.  Let's compare the two.  First, compute the opening using
% |imopen|.


se = strel('disk', 40);
Io = imopen(I, se);
%figure
%imshow(Io), title('Opening (Io)')

%%
% Next compute the opening-by-reconstruction using |imerode| and
% |imreconstruct|.


Ie = imerode(I, se);
Iobr = imreconstruct(Ie, I);
figure
imshow(Iobr), title('Opening-by-reconstruction (Iobr)')

%%
% Following the opening with a closing can remove the dark spots
% and stem marks.  Compare a regular morphological closing with a
% closing-by-reconstruction.  First try |imclose|:

Ioc = imclose(Io, se);
%figure
%imshow(Ioc), title('Opening-closing (Ioc)')

%%
% Now use |imdilate| followed by |imreconstruct|.  Notice you must 
% complement the image inputs and output of |imreconstruct|.

Iobrd = imdilate(Iobr, se);
Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
figure
imshow(Iobrcbr), title('Opening-closing by reconstruction (Iobrcbr)')

%%
% As you can see by comparing |Iobrcbr| with |Ioc|, reconstruction-based
% opening and closing are more 
% effective than standard opening and closing at removing small blemishes
% without affecting the overall shapes of the objects.  Calculate the
% regional maxima of |Iobrcbr| to obtain good foreground markers.


fgm = imregionalmax(Iobrcbr);
figure
imshow(fgm), title('Regional maxima of opening-closing by reconstruction (fgm)')

%%
% To help interpret the result, superimpose the foreground marker image
% on the original image. 

I2 = I;
I2(fgm) = 255;
%figure
%imshow(I2), title('Regional maxima superimposed on original image (I2)')

%%
% Notice that some of the mostly-occluded and shadowed objects are not
% marked, which means that these objects will not be segmented properly
% in the end result.  Also, the foreground markers in some objects go right
% up to the objects' edge.  That means you should clean the edges of the
% marker blobs and then shrink them a bit.  You can do this by a closing
% followed by an erosion.



se2 = strel(ones(5,5));
fgm2 = imclose(fgm, se2);
fgm3 = imerode(fgm2, se2);

%%
% This procedure tends to leave some stray isolated pixels that must be
% removed.  You can do this using |bwareaopen|, which removes all blobs
% that have fewer than a certain number of pixels. 


fgm4 = bwareaopen(fgm3, 20);
I3 = I;
I3(fgm4) = 255;
%figure
%imshow(I3)
%title('Modified regional maxima superimposed on original image (fgm4)')

%% Step 4: Compute Background Markers
% Now you need to mark the background.  In the cleaned-up image, |Iobrcbr|,
% the dark pixels belong to the background, so you could start with a
% thresholding operation.

bw = im2bw(Iobrcbr, graythresh(Iobrcbr)); %patient6.jpg at 0.82 , lambda 5.48 Iobrcbr
figure
imshow(bw), title('Thresholded opening-closing by reconstruction (bw)')

%%
% The background pixels are in black, but ideally we don't want the
% background markers to be too close to the edges of the objects we are
% trying to segment.  We'll "thin" the background by computing the
% "skeleton by influence zones", or SKIZ, of the foreground of |bw|.
% This can be done by computing the watershed transform of the distance
% transform of |bw|, and then looking for the watershed ridge lines (|DL ==
% 0|) of the result.
D = bwdist(bw);
DL = watershed(D);
bgm = DL == 0; %ridge lines
%figure
%imshow(bgm), title('Watershed ridge lines (bgm)')

%% Step 5: Compute the Watershed Transform of the Segmentation Function.
% The function |imimposemin| can be used to modify an image so that it has 
% regional minima only in certain desired locations.  Here you can use
% |imimposemin| to modify the gradient magnitude image so that its only
% regional minima occur at foreground and background marker pixels.
gradmag2 = imimposemin(gradmag, bgm | fgm4);

%%
% Finally we are ready to compute the watershed-based segmentation.
L = watershed(gradmag2);
Lrgb1 = label2rgb(L);
%figure,title('result of watershed transform'),imshow(Lrgb1);

%% Step 6: Visualize the Result
% One visualization technique is to superimpose the foreground
% markers, background markers, and segmented object boundaries on the
% original image.  You can use dilation as needed to make certain aspects,
% such as the object boundaries, more visible.  Object boundaries are
% located where |L == 0|.
I4 = I;
I4(imdilate(L == 0, ones(3, 3)) | bgm | bw) = 1;
figure
imshow(I4)
title('Markers and object boundaries superimposed on original image (I4)')

%%
% This visualization illustrates how the locations of the foreground and
% background markers affect the result.  In a couple of locations,
% partially occluded darker objects were merged with their brighter
% neighbor objects because the occluded objects did not have foreground
% markers.
%
% Another useful visualization technique is to display the label matrix
% as a color image.  Label matrices, such as those produced by
% |watershed| and |bwlabel|, can be converted to truecolor images for
% visualization purposes by using |label2rgb|.
%Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
%figure
%imshow(Lrgb)
%title('Colored watershed label matrix (Lrgb)')
%%
% You can use transparency to superimpose this pseudo-color label matrix on
% top of the original intensity image.
%figure
%imshow(I4)
%hold on
%himage = imshow(Lrgb);
%himage.AlphaData = 0.2;
%title('Lrgb superimposed transparently on original image')
%%
black=0;
white=0;
for i=1:(size(I4,1));
    for j=1:(size(I4,2));
        if I4(i,j)<=0.16
            black=black+1;
        else
            white=white+1;
        end
    end
end
threshold=17179;
black
white
if black>=threshold
    'Lung is normal'
else
    'Lung has Cancer'
end

%--------------------$

total1 = I4 - 0.3;
figure,imshow(total1),title('low brightness result');
total3 = medfilt2(total1); %preserves edges better than gaussian noise
%total3 = imgaussfilt(total1,1);
figure,imshow(total3),title('median filtered image');
total2 = im2bw(total3,0.65); %for patient1.jpg, 0.65 threshold
figure,imshow(total2)
total5 = bwareaopen(total2,150); %bwareaopen patient1.jpg = 150
figure,imshow(total5)
%separate nodule from lung wall
se1 = strel('square',3);  %for patient6.jpg,patient1.jpg 4 becomes 3 -se1 
totalerode = imerode(total5,se1);
figure,imshow(totalerode),title('final nodule');
left = totalerode(:, 1:end/2, :);
right = totalerode(:, end/2+1:end, :);
se2 = strel('square', 10);
right =imerode(right,se2);
figure,imshow(left),title('left lung');
figure,imshow(right),title('right lung');
%for left lung
Inoborder = imclearborder(left);
Imx = keepMaxObj(Inoborder);
figure,imshow(Imx);


%highlight nodule pixels
redAndBlueChannel = 255 * uint8(Imx);
% redChannel = 255* uint8(Imx); to generate red color
% blueChannel = 255*zeros(size(Imx),'uint8');
 greenChannel = 255 * zeros(size(Imx), 'uint8'); % Green Everywhere.
rgbImage1 = cat(3, redAndBlueChannel, greenChannel, redAndBlueChannel);
figure,imshow(rgbImage1);
figure
imshow(img_in)
hold on
himage = imshow(rgbImage1);
himage.AlphaData = 0.5;

stats = regionprops(Imx,'Area','Eccentricity','Perimeter'); 
area = stats.Area;
eccentricity = stats.Eccentricity;
eccentricity = im2double(eccentricity);
perimeter = stats.Perimeter;
perimeter = im2double(perimeter);
fprintf('Area is:  %d \n',area)
fprintf('Eccentricity is:  %d \n',eccentricity)
fprintf('Perimeter is: %d \n',perimeter)

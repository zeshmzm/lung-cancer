%% Preprocessing using gabor filter for image enhancement
lambda  = 6;
theta   = 0;
bw      = 3;
%sigma = lambda/pi*sqrt(log(2)/2)*(2^bw+1)/(2^bw-1);
psi     = [0 0];
gamma   = 1.95;
N       = 4;
path = 'F:\Local Disk D\Study\Final Year Project\Final Project\FYP 4th Feb9 watershed\Dataset\patient1.jpg';

%divide path into 3 parts for extension finding
[filepath,name,ext] = fileparts(path);

if strcmp(ext,'.dcm')
    img_in = dicomread(path);
else
      img_in = imread(path);
      img_in = imresize(img_in,[512 512]);
    
end

%img_in1 = imread(path);

%find if datatype is uint8 or not
if class(img_in) ~= 'uint8'
    img_in = im2uint8(img_in);
end 

%convert 3d into 2d
if size(img_in,3)==3
    img_in = rgb2gray(img_in);
end 
%find intensity values
i = img_in(:);  % image into 1 column
intensity = sum(i);
if intensity>20974604 && strcmp(ext,'.dcm')
    lambda = 5;
    gamma = 4;
    bw = 2.1;
    img_in = imadjust(img_in); %contrast stretching
elseif intensity>20974604
    lambda=5;
   
    
end

%img_in = im2double(img_i);
%img_in = dicomread('e.dcm');
%img_in(:,:,2:3) = [];


%gabor filter start
img_out = zeros(size(img_in,1), size(img_in,2), N);
for n=1:N
    gb = gabor_fn(bw,gamma,psi(1),lambda,theta)...
        + gabor_fn(bw,gamma,psi(2),lambda,theta);
    img_out(:,:,n) = imfilter(img_in, gb, 'symmetric');
    theta = theta + pi/4;
end
figure(1);
imshow(img_in);
title('input image');
img_out_disp = sum(abs(img_out).^2, 3).^0.5;

% right division
img_out_disp = img_out_disp./max(img_out_disp(:));
figure(2);
imshow(img_out_disp);
title('gabor output, L-2 super-imposed, normalized - imgoutdisp');

%% This is active contour segmentation
I = img_out_disp;
se = strel('disk', 20);
Ie = imerode(I, se);
figure(3),imshow(Ie),title('Eroded result')
Iobr = imreconstruct(Ie, I);
figure(4),imshow(Iobr),title('Reconstructed Eroded result')
Iobrd = imdilate(Iobr, se);
figure(5),imshow(Iobrd),title('Dilated result')
Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
figure(6),imshow(Iobrcbr),title('Reconstructed dilated result')
Iobrcbr = imcomplement(Iobrcbr);
figure(7),imshow(Iobrcbr),title(' Complemented reconstructred dilated result')
bw = im2bw(Iobrcbr, graythresh(Iobr));
figure
imshow(bw), title('Thresholded')
I=bw;
m = zeros(size(I,1),size(I,2));          %-- create initial mask
m(200:320,95:180) = 1; %1 black 0 white
m(186:321,348:410) = 1;
%I = imresize(I,.5);  %-- make image smaller 
%m = imresize(m,.5);  %     for fast computation
%subplot(2,2,1); imshow(I); title('Input Image');
%subplot(2,2,2); imshow(m); title('Initialization');
%subplot(2,2,3); title('Segmentation');
seg = region_seg(I, m, 1200); %-- Run segmentation
figure(8),imshow(seg), title('Global Region-Based Segmentation')
seg =im2uint8(seg);
I = seg;
se = strel('disk', 20);
Ie = imerode(I, se);
Iobr = imreconstruct(Ie, I);
Iobrd = imdilate(Iobr, se);
Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
bw = im2bw(Iobrcbr, graythresh(Iobr));
figure
imshow(bw), title('Thresholded')
 seg=bw;
 


%% Binarization for image classification
[M,N] = size(img_in);
total=ones(M,N);
white=0;
black=0;
for i=1:M;
    for j=1:N;
        if seg(i,j)==1
            total(i,j)=img_out_disp(i,j);
        else     
            total(i,j)=1;
        end
    end
end
figure(9),imshow(total),title('total image output')
for i=1:M;
    for j=1:N;
        if total(i,j)<=0.12
                black=black+1;%count black values
            else
                white=white+1;%count white values
        end
    end
end
figure,imshow(total);
threshold=17179
if black>=threshold
    ('normal lung')
else
    ('there is a possibility of nodules in the lungs')
end
total1 = total - 0.3;
figure,imshow(total1),title('low brightness result');
total3 = medfilt2(total1);
figure,imshow(total3),title('median filtered image');
total2 = im2bw(total3,0.55);
figure,imshow(total2)
total5 = bwareaopen(total2,500);
figure,imshow(total5)
%separate nodule from lung wall
se1 = strel('square',3);
totalerode = imerode(total5,se1);
figure,imshow(totalerode),title('final nodule');
left = totalerode(:, 1:end/2, :);
right = totalerode(:, end/2+1:end, :);
right =imerode(right,se1);
figure,imshow(left),title('left lung');
figure,imshow(right),title('right lung');

Inoborder = imclearborder(left);
Imx = keepMaxObj(Inoborder);
figure,imshow(Imx);

%highlight nodule pixels
redAndBlueChannel = 255 * uint8(Imx);
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


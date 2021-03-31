%% Preprocessing using gabor filter for image enhancement
lambda  = 7;
theta   = 0;
bw      = 3;
%sigma = lambda/pi*sqrt(log(2)/2)*(2^bw+1)/(2^bw-1);
psi     = [0 0];
gamma   = 2;
N       = 4;
path = 'a.bmp';
[filepath,name,ext] = fileparts(path);

if strcmp(ext,'.dcm')
    img_in = dicomread(path);
else
    img_in = imread(path);
    %img_in = imresize(img_in,[512 512]);
    
end

%img_in1 = imread(path);

if class(img_in) ~= 'uint8'
    img_in = im2uint8(img_in);
end 
if size(img_in,3)==3
    img_in = rgb2gray(img_in);
end 
i = img_in(:);
intensity = sum(i);

if intensity>9524455 && strcmp(ext,'.dcm')
    lambda = 5;
    gamma = 4;
    bw = 2.1;
    img_in = imadjust(img_in);
elseif intensity>9524455
    lambda=5;
end

%img_in = im2double(img_i);
%img_in = dicomread('e.dcm');
%img_in(:,:,2:3) = [];
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
figure(2);
img_out_disp = sum(abs(img_out).^2, 3).^0.5;
img_out_disp = img_out_disp./max(img_out_disp(:));
imshow(img_out_disp);
title('gabor output, L-2 super-imposed, normalized');

%% This is marker controlled watershed using masking
I = img_out_disp;
se = strel('disk', 20);
Ie = imerode(I, se);
Iobr = imreconstruct(Ie, I);
Iobrd = imdilate(Iobr, se);
Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
bw = im2bw(Iobrcbr, graythresh(Iobr));
figure
imshow(bw), title('Thresholded')
I=bw;
m = zeros(size(I,1),size(I,2));          %-- create initial mask
m(200:320,95:180) = 1;
m(186:321,348:410) = 1;
%I = imresize(I,.5);  %-- make image smaller 
%m = imresize(m,.5);  %     for fast computation
%subplot(2,2,1); imshow(I); title('Input Image');
%subplot(2,2,2); imshow(m); title('Initialization');
%subplot(2,2,3); title('Segmentation');
seg = region_seg(I, m, 1200); %-- Run segmentation
figure
imshow(seg); title('Global Region-Based Segmentation')
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
hasil=ones(M,N);
white=0;
black=0;
for i=1:M;
    for j=1:N;
        if seg(i,j)==1
            hasil(i,j)=img_out_disp(i,j);
            
        else
            hasil(i,j)=1;
        end
    end
end
for i=1:M;
    for j=1:M;
        if hasil(i,j)<=0.12
                black=black+1;
            else
                white=white+1;
        end
    end
end
figure,imshow(hasil);
threshold=17179
if black>=threshold
    ('normal lung')
else
    ('there is a possibility of nodules in the lungs')
end
hasil1 = hasil - 0.3;
figure,imshow(hasil1),title('low brightness result');
hasil3 = medfilt2(hasil1);
figure,imshow(hasil3),title('median filtered image');
hasil2 = im2bw(hasil3,0.55);
figure,imshow(hasil2)
hasil5 = bwareaopen(hasil2,500);
figure,imshow(hasil5)
se1 = strel('square',3);
hasilerode = imerode(hasil5,se1);
figure,imshow(hasilerode),title('final nodule');
left = hasilerode(:, 1:end/2, :);
right = hasilerode(:, end/2+1:end, :);
right =imerode(right,se1);
figure,imshow(left),title('left lung');
figure,imshow(right),title('right lung');

stats = regionprops(left,'Area','Eccentricity','Perimeter'); 
area = stats(2).Area;
eccentricity = stats(2).Eccentricity;
eccentricity = im2double(eccentricity);
perimeter = stats(2).Perimeter;
perimeter = im2double(perimeter);
fprintf('Area is:  %d \n',area)
fprintf('Eccentricity is:  %d \n',eccentricity)
fprintf('Perimeter is: %d \n',perimeter)


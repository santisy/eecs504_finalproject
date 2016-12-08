function siftRaw = siftDes(column,row,scale,image)
% This function generates a sift-histogram descriptor
% column and row are the coordinate of the feature point
% scale is the sigma for the DoG operator
% image is the original image
% sift part and histogram part will be normalized seperately
% the output descriptor has the following sturcture: column 1:2 coordinate,
% column 2:x: sift, column x:end historgram

% now start the sift part
% firstly get the image patch corresponding to patch
% filters

NUM_SUB_REGION = 16;

grayImage = rgb2gray(image);
[m,n] = size(grayImage);
keypointOrient = rad2deg(atan2((double(grayImage(row,column+1))-double(grayImage(row,column-1))),double(grayImage(row+1,column))-double(grayImage(row-1,column))));
% now rotate the image and fetch the real patch for calculation
T1 = maketform('affine',[1 0 0; 0 1 0; -column -row 1]);
R1 = maketform('affine',[cosd(-keypointOrient) sind(-keypointOrient) 0; -sind(-keypointOrient) cosd(-keypointOrient) 0; 0 0 1]);
T2 = maketform('affine',[1 0 0; 0 1 0; column row 1]);
transMat = maketform('composite', T2, R1, T1);
rotatedImage = imtransform(grayImage, transMat, 'XData', [1 m], 'YData', [1 n]);
patch = rotatedImage(floor(row-3*scale):floor(row+3*scale),floor(column-3*scale):floor(column+3*scale));
patch = imresize(patch,[16,16]); % get the rescaled path

% now calculate gradient direction and magnitude
dxFilter = [1,-1];
dyFilter = [1,-1]';
dxMat = double(imfilter(patch, dxFilter,'same'));
dyMat = double(imfilter(patch,dyFilter,'same'));
magnMat = sqrt(dxMat.*dxMat + dyMat.*dyMat);
orientMat = atan2(dyMat, dxMat);
% then do the gaussian weighting
magnMat = magnMat.*gaussianKernel(16);
siftRaw = zeros(16,8)
% then calculate the histogram
for i = 1:4:16
    for j = 1:4:16
        siftRaw(8*(i-1)+j) = hist(magnMat(i:i+3,j:i+3),orientMat(i:i+3,j:j+3));
    end
end
siftRaw = reshape(siftRaw,[128,1])
siftRaw = siftRaw./norm(siftRaw)

function h = gaussianKernel(size)
sigma = size./2;
ind = -floor(size/2) : floor(size/2);
[X Y] = meshgrid(ind, ind);
h = exp(-(X.^2 + Y.^2) / (2*sigma*sigma));
h = imresize(h,[16,16])
h = h / sum(h(:));

function result = hist(subMagn,subOrient) % eight bins
% orientation is the orientation ranging from -pi to pi, in rad
% subMagn is the magnitude
% I was wondering how to 
bins = linspace(-pi,pi,9);
% this will not be fast
offset = 0.5*(bins(4)-bins(3));
bins = bins+offset;
bins(end) = bins(end) - offset;
[m,n] = size(subOrient);
hist = zeros(8,1);
for i = 1:m
    for j = 1:n
        temp = bins - subOrient(i,j);
        for k = 1:7
            if temp(k)<=0 && temp(k+1)>=0
                hist(k) = hist(k) + subMagn(i,j);
                break
            end                    
        end
        hist(8) = hist(8) + subMagn(i,j);   
    end
end


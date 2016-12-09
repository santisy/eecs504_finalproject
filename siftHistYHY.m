function descriptor = siftHistYHY(row,column,scale,image)
% This function generates a sift-histogram descriptor
% Modified to be compatible with Yuhao Yu's script
% column and row are the coordinate of the feature point
% scale is the sigma for the DoG operator
% image is the original image
% sift part and histogram part will be normalized seperately
% the output descriptor has the following sturcture: column 1:30, histogram part
% column 31:158 sift part
% histogram and sift part are normalized seperately

% now start the sift part
%disp('debug')
image = imfilter(image, fspecial('gaussian'),'same');
grayImage = rgb2gray(image);
[m,n] = size(grayImage);
% now calculate the dominant orientation
dxFilter = [1,0,-1];
dyFilter = [1,0,-1]';

prePatch = grayImage(floor(row-3*scale):floor(row+3*scale),floor(column-3*scale):floor(column+3*scale),:);
preDxMat = double(imfilter(prePatch, dxFilter,'symmetric'));
preDyMat = double(imfilter(prePatch,dyFilter,'symmetric'));
preMagnMat = sqrt(preDxMat.*preDxMat + preDyMat.*preDyMat);
[m1,n1] = size(preDxMat);
% preDxMat = gaussianKernel(m1).*preDxMat;
% preDyMat = gaussianKernel(m1).*preDyMat;
% Ixy = [preDxMat(:)';preDyMat(:)'];
% sTensor = Ixy*Ixy'; 
% [U,lbd] = eig(sTensor);
% u = U(:,1);
preOrientMat = atan2(preDyMat, preDxMat);
[keypointOrient,~] = gradHist(gaussianKernel(m1).*preMagnMat,preOrientMat,36);
%keypointOrient = rad2deg(keypointOrient);



%keypointOrient = rad2deg(atan2((double(grayImage(row+1,column))-double(grayImage(row-1,column))),double(grayImage(row,column+1))-double(grayImage(row,column-1))));
% now rotate the image and fetch the real patch for calculation
T1 = maketform('affine',[1 0 0; 0 1 0; -column -row 1]);
R1 = maketform('affine',[cos(-keypointOrient) sin(-keypointOrient) 0; -sin(-keypointOrient) cos(-keypointOrient) 0; 0 0 1]);
T2 = maketform('affine',[1 0 0; 0 1 0; column row 1]);
transMat = maketform('composite', T2, R1, T1);
rotatedImage = imtransform(image, transMat, 'XData', [1 m], 'YData', [1 n]);
patch = rotatedImage(floor(row-3*scale):floor(row+3*scale),floor(column-3*scale):floor(column+3*scale),:);
patch = imresize(patch,[16,16]); % get the rescaled path

histVec = zeros(30,1);
% now calculate the hist
centers = 1/10;       % bin offset
colorBins = (centers/2:centers:1)*128;   % bin centers
for i = 1:3 %iteration through channels
    temp = patch(:,:,i);
    histTemp = hist(temp(:),colorBins);
    histVec((i-1)*10+1:i*10) = histTemp;
end
histVec = histVec./norm(histVec);

% now calculate gradient direction and magnitude
patch = rgb2gray(patch);
dxMat = double(imfilter(patch, dxFilter,'symmetric'));
dyMat = double(imfilter(patch,dyFilter,'symmetric'));
magnMat = sqrt(dxMat.*dxMat + dyMat.*dyMat);
orientMat = atan2(dyMat, dxMat);
% then do the gaussian weighting
magnMat = magnMat.*gaussianKernel(16);
siftRaw = zeros(16,8);
% then calculate the histogram
count = 0;
for i = 1:4:16
    for j = 1:4:16
        count = count + 1;
        [~,siftRaw(count,:)] = gradHist(magnMat(i:i+3,j:j+3),orientMat(i:i+3,j:j+3),8);
    end
end
siftRaw = reshape(siftRaw,[128,1]);
siftRaw = siftRaw./norm(siftRaw);

descriptor = zeros(30+128,1);
descriptor(1:30) = histVec;
descriptor(31:end) = siftRaw;

%now conbine all the elements

function h = gaussianKernel(size)
sigma = size./2;
ind = -floor(size/2) : floor(size/2);
[X Y] = meshgrid(ind, ind);
h = exp(-(X.^2 + Y.^2) / (2*sigma*sigma));
h = imresize(h,[size,size]);
h = h / sum(h(:));

function [maxOrient,hist] = gradHist(subMagn,subOrient,numBin) % eight bins
% orientation is the orientation ranging from -pi to pi, in rad
% subMagn is the magnitude
bins = linspace(-pi,pi,numBin+1);% bin wall
% this will not be fast
offset = 0.5*(bins(4)-bins(3));% 
[m,n] = size(subOrient);
hist = zeros(numBin,1);
binCenter = bins+offset;
for i = 1:m
    for j = 1:n
        temp = bins - subOrient(i,j);
        flag = 0;
        for k = 1:(numBin-1)
            if temp(k)<=0 && temp(k+1)>=0
                hist(k+1) = hist(k+1) + subMagn(i,j);
                flag = 1;
                break
            end                    
        end
        if flag == 0
            hist(1) = hist(1) + subMagn(i,j);
        end
    end
end
[~,maxIdx] = max(hist);
maxOrient = binCenter(maxIdx);
if maxOrient >= pi
    maxOrient = maxOrient - 2*pi;
end


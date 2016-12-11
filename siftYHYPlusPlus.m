function descriptor = siftHistYHYPlus(row,column,scale,keypointOrient,image)

% This function generates a sift-histogram descriptor
% The orientation is given by the "orientation" term
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
image = imfilter(image, fspecial('gaussian',[8,8],scale),'same');
grayImage = rgb2gray(image);
[m,n] = size(grayImage);
% % now calculate the dominant orientation
dxFilter = [1,0,-1];
dyFilter = [1,0,-1]';
% 
% prePatch = grayImage(floor(row-3*scale):floor(row+3*scale),floor(column-3*scale):floor(column+3*scale),:);
% preDxMat = double(imfilter(prePatch, dxFilter,'symmetric'));
% preDyMat = double(imfilter(prePatch,dyFilter,'symmetric'));
% preMagnMat = sqrt(preDxMat.*preDxMat + preDyMat.*preDyMat);
% [m1,n1] = size(preDxMat);
% preDxMat = gaussianKernel(m1).*preDxMat;
% preDyMat = gaussianKernel(m1).*preDyMat;
% Ixy = [preDxMat(:)';preDyMat(:)'];
% sTensor = Ixy*Ixy'; 
% [U,lbd] = eig(sTensor);
% u = U(:,1);
% preOrientMat = atan2(preDyMat, preDxMat);
% [keypointOrient,~] = gradHist(gaussianKernel(m1).*preMagnMat,preOrientMat,36);
%keypointOrient = rad2deg(keypointOrient);

%keypointOrient = rad2deg(atan2((double(grayImage(row+1,column))-double(grayImage(row-1,column))),double(grayImage(row,column+1))-double(grayImage(row,column-1))));
% now rotate the image and fetch the real patch for calculation
T1 = maketform('affine',[1 0 0; 0 1 0; -column -row 1]);
R1 = maketform('affine',[cos(-keypointOrient) sin(-keypointOrient) 0; -sin(-keypointOrient) cos(-keypointOrient) 0; 0 0 1]);
T2 = maketform('affine',[1 0 0; 0 1 0; column row 1]);
transMat = maketform('composite', T2, R1, T1);
rotatedImage = imtransform(image, transMat, 'XData', [1 m], 'YData', [1 n]);
% patch = rotatedImage(floor(row-3*scale):floor(row+3*scale),floor(column-3*scale):floor(column+3*scale),:);
% patch = imresize(patch,[16,16]); % get the rescaled path
 patch = rotatedImage(floor(row-0.5-8):floor(row - 0.5+8),floor(column-0.5-8):floor(column - 0.5+8),:);
 patch = imresize(patch,[16,16]);
%  close all;
%  imagesc(image);
%  hold on
%  plot(column,row,'r*');
%  figure
%  imagesc(rotatedImage);
%  hold on
%  plot(column,row,'r*');
 
 
% histVec = zeros(30,1);
% now calculate the hist
% centers = 1/10;       % bin offset
% colorBins = (centers/2:centers:1)*256;   % bin centers
% for i = 1:3 %iteration through channels
%     temp = patch(:,:,i);
%     histTemp = hist(temp(:),colorBins);
%     histVec((i-1)*10+1:i*10) = histTemp;
% end
% histVec = histVec./norm(histVec);


% now calculate gradient direction and magnitude
patch = im2double(rgb2gray(patch));
dxMat = (imfilter(patch, dxFilter,'symmetric'));
dyMat = (imfilter(patch, dyFilter,'symmetric'));
magnMat = sqrt(dxMat.*dxMat + dyMat.*dyMat);
% figure
% imagesc(patch);
% figure
% imagesc(dxMat);
% figure
% imagesc(dyMat);
% figure
% imagesc(magnMat);

magnMat = magnMat.*gaussianKernel(16,0.5*16);
% figure
% imagesc(magnMat);
orientMat = atan2(dyMat, dxMat);
siftRaw = grad_hist(magnMat,orientMat,8,scale);

siftRaw = siftRaw./norm(siftRaw);
idx = find(siftRaw>=0.2);
siftRaw(idx) = 0.2;
siftRaw = siftRaw./norm(siftRaw);

%descriptor = zeros(30+128,1);
descriptor = siftRaw;
%descriptor(1:30) = histVec;
%descriptor(31:end) = siftRaw;

%now conbine all the elements

function h = gaussianKernel(size,sigma)
ind = -floor(size/2) : floor(size/2);
[X Y] = meshgrid(ind, ind);
h = exp(-(X.^2 + Y.^2) / (2*sigma*sigma));
h = imresize(h,[size,size]);
h = h / sum(h(:));

function sub_hist = sub_Hist(subMagn,subOrient,numBin) % eight bins
% orientation is the orientation ranging from -pi to pi, in rad
% subMagn is the magnitude
subMagn = subMagn(:);
subOrient = subOrient(:);
bins = linspace(-pi,pi,numBin+1);% bin wall
hist_count = zeros(1,numBin);
num_of_pix = length(subMagn);
sub_temp = subOrient*ones(1,numBin+1)- ones(num_of_pix,1)*bins;
sub_temp(sub_temp<=0) = 0;
sub_temp(sub_temp>0) = 1;
index_bin = sum(sub_temp,2);
for k = 1:numBin
    hist_count(k) = sum(subMagn(index_bin==k));
end
sub_hist = hist_count;

function descriptor = grad_hist(magn_patch,orient_patch,num_bin,scale)
% assuming the patch is 16*16
% And the image has been rotated
descriptor = [];
for i = 1:4:16
    for j = 1:4:16
        subMagn = magn_patch(i:i+3,j:j+3);
        subOrient = orient_patch(i:i+3,j:j+3);
        subHist = sub_Hist(subMagn,subOrient,num_bin);
        descriptor = [descriptor,subHist];
    end
end

        
        



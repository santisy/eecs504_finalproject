% This code is based on VLFeat open source lib
% Before running this code, make sure VLFeat Lib is installed
% run('VLFEATROOT/toolbox/vl_setup'); VLFEATROOT denotes the VLFEAT
% directory, e.g. '~/Downloads/VLFeat'
% The code does three things listed as below
% 1. Feature points detector : Harris Laplace
% 2. Feature descriptor : SIFT and patch color histogram
close all
clear
load 'rect_info.mat'
% read two image
image1id = 2;
image2id = 3;
IM1 = imread(['00',num2str(image1id),'.png']);
I1 = single(rgb2gray(IM1));

IM2 = imread(['00',num2str(image2id),'.png']);
I2 = single(rgb2gray(IM2));

% smooth

H = fspecial('gaussian');
I1 = conv2(I1,H,'same');
I2 = conv2(I2,H,'same');

%--------------------------------------------------------------------------
%% Feature detection
%--------------------------------------------------------------------------
% Perform Harris Laplace detector
% vl_covdet returns feature points, SIFT descriptors, and 
% feature measurements 
[F1,D1,INFO1] = vl_covdet(I1,'Method','HarrisLaplace');
[F2,D2,INFO2] = vl_covdet(I2,'Method','HarrisLaplace');
% cornerness measurement 
peaktol = 10;
peak_mask1 = abs(INFO1.peakScores) >= peaktol;
peak_mask2 = abs(INFO2.peakScores) >= peaktol;
% scale space peak measurement
scaletol = 5;
scale_mask1 = abs(INFO1.laplacianScaleScore) >= scaletol;
scale_mask2 = abs(INFO2.laplacianScaleScore) >= scaletol;

mask1 = peak_mask1 & scale_mask1;
mask2 = peak_mask2 & scale_mask2;
% filter out weak features
F1 = F1(:,mask1);
D1 = D1(:,mask1);
F2 = F2(:,mask2);
D2 = D2(:,mask2);
%--------------------------------------------------------------------------
%% Group points in rectanglar
%--------------------------------------------------------------------------
table1 = group_points(image1id,rect_record,F1);
F1 = F1(:,table1(:,1));
D1 = D1(:,table1(:,1));
point_record1 = table1(:,2); % record rectanglar index

table2 = group_points(image2id,rect_record,F2);
F2 = F2(:,table2(:,1));
D2 = D2(:,table2(:,1));
point_record2 = table2(:,2);
%--------------------------------------------------------------------------
%% Build color histogram
%--------------------------------------------------------------------------
% extract patches
R = 5; % set standard patch size as 21-by-21
R1 = floor(F1(3,:)+0.5);
xc = F1(1,:);
yc = F1(2,:);
n1 = size(F1,2); % number of keypoints in im1
P1 = zeros((2*R+1)^2,n1,3);
for i = 1:n1
    % Assume we would never cross the boundary
    xlt = floor(xc(i) - R1(i) + 0.5);
    ylt = floor(yc(i) - R1(i) + 0.5);
    pat = IM1(ylt:ylt+2*R1(i),xlt:xlt+2*R1(i),:);
    % normalize patch size to standard size
    pat = imresize(pat,[2*R+1,2*R+1]); 
    P1(:,i,:) = reshape(pat,[(2*R+1)^2,1,3]);
end
R2 = floor(F2(3,:)+0.5);
xc = F2(1,:);
yc = F2(2,:);
n2 = size(F2,2); % number of keypoints in im2
P2 = zeros((2*R+1)^2,length(R2),3);
for i = 1:length(R2)
    % Assume we would never cross the boundary
    xlt = floor(xc(i) - R2(i) + 0.5);
    ylt = floor(yc(i) - R2(i) + 0.5);
    pat = IM2(ylt:ylt+2*R2(i),xlt:xlt+2*R2(i),:);
    pat = imresize(pat,[2*R+1,2*R+1]);
    P2(:,i,:) = reshape(pat,[(2*R+1)^2,1,3]);
end

% build color histgram
b = 20;       % 128 bins
c = 1/b;       % bin offset
x = c/2:c:1;   % bin centers

DC1 = zeros(3*b,n1); % color histogram descriptor
for i = 1:n1
    v = [];
    for j = 1:3 
        h = hist(single(P1(:,i,j))/255,x)';
        v = [v;h/sum(h)]; % noramlize histgram vector
    end
    DC1(:,i) = v;
end

DC2 = zeros(3*b,n2); % color histogram descriptor
for i = 1:n2
    v = [];
    for j = 1:3 
        h = hist(single(P2(:,i,j))/255,x)';
        v = [v;h/sum(h)]; % noramlize histgram vector
    end
    DC2(:,i) = v;
end

DC1 = single(DC1);
DC2 = single(DC2);
%--------------------------------------------------------------------------
%% Extract histogram in rectangular
bin = 20;
rectlist = rect_record(rect_record(:,9) == image1id,:);
V1 = zeros(3*bin,size(rectlist,1));
for i = 1:size(rectlist,1)
    v = rect_hist(IM1,rectlist(i,:),bin);
    V1(:,i) = v/norm(v);
end

rectlist = rect_record(rect_record(:,9) == image2id,:);
V2 = zeros(3*bin,size(rectlist,1));
for i = 1:size(rectlist,1)
    v = rect_hist(IM2,rectlist(i,:),bin);
    V2(:,i) = v/norm(v);
end
%--------------------------------------------------------------------------
%% Match
%--------------------------------------------------------------------------
% default threshold = 1.5
% distance between matched points <= thresh * dist between the point and
% anyother points
thresh = 1.5;
rect_num1 = length(unique(point_record1));
rect_num2 = length(unique(point_record2));
cost = zeros(rect_num1,rect_num2);
De1 = [D1;1.5*DC1];
De2 = [D2;1.5*DC2];
for i = 1:rect_num1
    idx1 = (point_record1 == i);
    for j = 1:rect_num2     
        idx2 = (point_record2 == j);
        %{
        [~,score1] = fastmatch(D1(:,idx1),D2(:,idx2));
        [~,score2] = fastmatch(2*DC1(:,idx1),2*DC2(:,idx2));
        score1 = sort(score1,'ascend');
        score2 = sort(score2,'ascend');
        %}
        
        [~,score] = fastmatch(De1(:,idx1),De2(:,idx2),1.1);
        %[~,score] = vl_ubcmatch(De1(:,idx1),De2(:,idx2),thresh);
        score = sort(score,'ascend');
        if (length(score) < 5)
            cost(i,j) = inf;
        else
            cost(i,j) = mean(score(1:5));
        end
        %cost(i,j) = cost(i,j) + norm(V1(:,i)-V2(:,j))/2;
    end
end
%--------------------------------------------------------------------------
%% Show match
%--------------------------------------------------------------------------
% 
[~,idx] = min(cost(:));
[i,j] = ind2sub(size(cost),idx);
rectlist1 = rect_record(rect_record(:,9) == image1id,:);
rect_info1 = rectlist1(i,:);
rectlist2 = rect_record(rect_record(:,9) == image2id,:);
rect_info2 = rectlist2(j,:);
imagesc(IM1), hold on;
draw_rect(rect_info1);
figure;
imagesc(IM2), hold on;
draw_rect(rect_info2);
%{
reorder 
%[scores,idx] = sort(scores,'descend');
%matches = matches(:,idx);
% matching points will be indicated by the same number
nm = size(matches,2);
figure;
imagesc(IM1), hold on;
vl_plotframe(F1(:,matches(1,1:nm)));
for i = 1:nm
    x = F1(1,matches(1,i));
    y = F1(2,matches(1,i));
    text(x,y,num2str(i));
end
hold off

figure;
imagesc(IM2), hold on;
vl_plotframe(F2(:,matches(2,1:nm)));
for i = 1:nm
    x = F2(1,matches(2,i));
    y = F2(2,matches(2,i));
    text(x,y,num2str(i));
end
%}

% test group pixels 
close all
load 'rect_info'
image1id = 3;
image2id = 5;
IM1 = imread(['00',num2str(image1id),'.png']);
I1 = single(rgb2gray(IM1));

IM2 = imread(['00',num2str(image2id),'.png']);
I2 = single(rgb2gray(IM2));

% feature measurements 
%[F1,D1,INFO1] = vl_covdet(I1,'Method','HarrisLaplace');
%[F2,D2,INFO2] = vl_covdet(I2,'Method','HarrisLaplace');
pt1 = load(['points00',num2str(image1id),'.mat']);
pt2 = load(['points00',num2str(image2id),'.mat']);
F1 = pt1.pt';
F1 = [F1;zeros(1,size(F1,2))];
F1(1:2,:) = flipud(F1(1:2,:));
F2 = pt2.pt';
F2 = [F2;zeros(1,size(F2,2))];
F2(1:2,:) = flipud(F2(1:2,:));

% build SIFT descriptor
peakthresh = 2;
[F1,D1] = vl_sift(I1,'Frames',F1,'Orientations','PeakThresh',peakthresh);
[F2,D2] = vl_sift(I2,'Frames',F2,'Orientations','PeakThresh',peakthresh);

% group 
table1 = group_points(image1id,rect_record,F1);
table2 = group_points(image2id,rect_record,F2);

% show
subplot 121
imagesc(IM1), hold on;
vl_plotframe(F1);
subplot 122
show_pointinrect(table1,F1,image1id,rect_record);

figure;
subplot 121
imagesc(IM2), hold on;
vl_plotframe(F2);
subplot 122
show_pointinrect(table2,F2,image2id,rect_record);




clc,clear
close all
load('002_circle.mat');
bw_im = zeros(400,400,5);
for i = 1:5
    bw_im(:,:,i) = imopen(imbinarize(rgb2gray(double(imread(['00', num2str(i),'.png']))/255),'adaptive','ForegroundPolarity',...
                    'dark','Sensitivity',0.5),strel('disk',5));
end
bw_im(~repmat(B_circle, [1,1,5])) = 1;
bw_im = ~bw_im;
%imshow(bw_im(:,:,1)) 
rect_record = [];% this will be a n x 9 rectangular record, length direction and
                 % width direction are respectively 1 x 2 variable
                 % [x, y, length/2, length_direction, width/2, width_direction,pic#]
%label_bw = zeros(size(bw_im));
rect_count = 1;
for i = 1:5
    label_bw = bwlabel(bw_im(:,:,i));
    unique_label = unique(label_bw);
    for j = 1:length(unique_label)
        if sum(sum(label_bw==unique_label(j)))<15
            label_bw(label_bw==unique_label(j)) = 0;
        end
    end
    unique_label = unique(label_bw);
    unique_label(unique_label==0) = [];
    % if unique_label number greater than 8, we want to delete the rest few
    % because the one may not be sensitive to detect. Here is the upper
    % side of the ice cube shape
    if length(unique_label)>8
        num_of_pix = zeros(1, length(unique_label));
        for j = 1:length(unique_label)
            num_of_pix(j) = sum(sum(label_bw==unique_label(j)));
        end
        [~,index_descend] = sort(num_of_pix,'descend');
        unique_label = unique_label(index_descend(1:8));
    end
    rect_tempRecord = zeros(length(unique_label), 8);
    for j = 1:length(unique_label)
        [lx, ly] = ind2sub([400,400], find(label_bw==unique_label(j)));
        point_l = [lx, ly];
        rect_tempRecord(j,1:2) = mean(point_l);
        point_l_til = point_l - ones(length(point_l),1)*mean(point_l);
        [Q,D] = eig(point_l_til'*point_l_til);
        rect_tempRecord(j,7:8) = Q(:,1)';
        rect_tempRecord(j,4:5) = Q(:,2)';
        point_l_Q = Q'*point_l_til';
        [~,index_length]= max(point_l_Q(2,:));
        [~,index_width] = max(point_l_Q(1,:));
        rect_tempRecord(j,3) = norm(point_l(index_length,:)-mean(point_l),2);
        rect_tempRecord(j,6) = norm(point_l(index_width,:)-mean(point_l),2);
    end
    rect_record(rect_count:(rect_count+length(unique_label)-1),1:8) = rect_tempRecord;
    rect_record(rect_count:(rect_count+length(unique_label)-1), 9) = ones(length(unique_label),1)*i;
    rect_count = rect_count + length(unique_label);
end

for display_key = 1:5
	I = imread(['00',num2str(display_key),'.png']);
	figure()
	imshow(I);
	hold on
	draw_rect(rect_record(rect_record(:,end)==display_key,:));
	hold off
	figure()
	imshow(bw_im(:,:,display_key))
end




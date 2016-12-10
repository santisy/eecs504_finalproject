function [pointLoc,scales,cornerScore,scaleScore,orientation] = harrislaplacePlusPlus(img,corner_threshold,scale_thresold)
    % modified version dec 09
    % modified version dec 08
    % orientation term will define the dominant orientation of the image
    % patch for further SIFT descriptor
    % A single point with 2 dominant directions will be recognized as 2
    % points
    % the score for cornerness and scaleness are also returned by the 
    % harris laplace
    % the ratio of points to be reserved based on cornerness is controlled
    % by conernessRatio (0~1), as a threshold, the higher this value is,
    % the less points would be preserved
    % based on this person's work, the original file name is
    % kp_harrislaplace
    % Author :: Vincent Garcia
    % Date   :: 05/12/2007
    
    % this file can provide points to describe in SIFT.
    
    % this file take in gray picture and return points, 
    % every row of the points variable is a detected key point, the first
	% two are the rows and columns(all in matrix way) NOTE: if you 
	% want to plot, you have to follow picture convention, where you have to
	% converse your coordinates
    % What I have done to the file:
    % 1)modify some coefficient of the algorithm to the adapt to the
    % specific situation
    % 2)the 002_circle.mat is a precalculated circle range
    % 3)the finally detected key points, actually, many are very close, so
    % I finally, combine close points, the scale of them are the minimum
    % one among them(see last few lines of code in this file)
    
    if corner_threshold == []
        corner_threshold = 0.5 % by default
    end
    if scale_thresold == []
        scale_thresold = 0.5 % by default
    end
    % IMAGE PARAMETERS
    load('002_circle.mat');
    img         = double(img(:,:,1));
    img_height  = size(img,1);
    img_width   = size(img,2);

    % SCALE PARAMETERS
    sigma_begin = 0.5;
    sigma_step  = 1.2;
    sigma_nb    = 12;
    sigma_array = (sigma_step.^(0:sigma_nb-1))*sigma_begin;
    k = 0.06;
    max_ratio = corner_threshold;
    % PART 1 : HARRIS
    harris_pts = zeros(0,6);% 4th column: corner score 5th column: laplace response score 6th column: orientation
   
    for i=1:sigma_nb

        % scale (standard deviation)
        s_I = sigma_array(i);   % integration scale
        s_D = 0.7*s_I;          % derivative scale %0.7

        % derivative mask
        x  = -round(3*s_D):round(3*s_D);
        dx = x .* exp(-x.*x/(2*s_D*s_D)) ./ (s_D*s_D*s_D*sqrt(2*pi));
        dy = dx';

        % image derivatives
        Ix = conv2(img, dx, 'same');
        Iy = conv2(img, dy, 'same');

        % auto-correlation matrix
        g   = fspecial('gaussian',max(1,fix(6*s_I+1)), s_I);
        Ix2 = conv2(Ix.^2, g,  'same');
        Iy2 = conv2(Iy.^2, g,  'same');
        Ixy = conv2(Ix.*Iy, g, 'same');

        % interest point response
        %cim = (Ix2.*Iy2 - Ixy.^2)./(Ix2 + Iy2 + eps);				% Alison Noble measure.
        cim = (Ix2.*Iy2 - Ixy.^2) - k*(Ix2 + Iy2).^2;	% Original Harris measure.
        cim(~B_circle) = 0;
        % find local maxima on neighborgood
        [~,~,max_local] = findLocalMaximum(cim,3*s_I);%3*s_I
        % set threshold 1% of the maximum value
        t = max_ratio*max(max_local(:));

        % find local maxima greater than threshold
        [l,c] = find(max_local>=t);

        % build interest points
        n = size(l,1);
        harris_pts(end+1:end+n,:) = [l,c,repmat(i,[n,1]),max_local(sub2ind([400,400],l,c)),repmat(0,[n,1]),repmat(0,[n,1])];
    end


    % PART 2 : LAPLACE
    % compute scale-normalized laplacian operator
    laplace_snlo = zeros(img_height,img_width,sigma_nb);
    for i=1:sigma_nb
        s_L = sigma_array(i);   % scale
        laplace_snlo(:,:,i) = s_L*s_L*imfilter(img,fspecial('log', floor(6*s_L+1), s_L),'replicate');
    end
    % verify for each of the initial points whether the LoG attains a maximum at the scale of the point
    n   = size(harris_pts,1);
    cpt = 0;
    points = zeros(n,6);
    for i=1:n
        l = harris_pts(i,1);
        c = harris_pts(i,2);
        s = harris_pts(i,3);
        val = laplace_snlo(l,c,s);
        if s>1 && s<sigma_nb
            if val>laplace_snlo(l,c,s-1) && val>laplace_snlo(l,c,s+1)
                cpt = cpt+1;
                points(cpt,:) = harris_pts(i,:);
                points(cpt,5) = val;
            end
        elseif s==1
            if val>laplace_snlo(l,c,2)
                cpt = cpt+1;
                points(cpt,:) = harris_pts(i,:);
                points(cpt,5) = val;
            end
        elseif s==sigma_nb
            if val>laplace_snlo(l,c,s-1)
                cpt = cpt+1;
                points(cpt,:) = harris_pts(i,:);
                points(cpt,5) = val;
            end
        end
    end

    points(cpt+1:end,:) = [];
    
    % Filters the laplacian response with threshold scaleRatio
    points_temp = points;
    max_response = max(points(:,5));
    threshold = max_response*scale_thresold;
    points = points_temp(find(abs(points_temp(:,5))>=threshold),:);
    
    % SET SCALE TO 3*SIGMA FOR DISPLAY
    points(:,3) = sigma_array(points(:,3));
    %% finally start to eliminate too close points
    m = size(points, 1);
    record_list = zeros(m,1);
    label = 1;
    new_points = [];
    while any(record_list==0)
        tdp = find(record_list==0,1);
        dist_list = abs(points(:,1:2) - repmat(points(tdp,1:2), [m,1]));
        record_list( dist_list(:,1)<=2&dist_list(:,2)<=2) = label;
        new_points(label,1:2) = points(tdp,1:2);
        new_points(label,4:6) = points(tdp,4:6);
        new_points(label,3) = min(points(dist_list(:,1)<=2&dist_list(:,2)<=2, 3));
        
        label = label + 1;   
    end
    % new_points;
    % STEP 4: ASSIGN DOMINANT DIRECTIONS
    points = zeros(0,6);
    magn_map = sqrt(Ix.*Ix + Iy.*Iy);
    orient_map = atan2(Iy, Ix);
    % then iteration through all the points
    [count,~] = size(new_points);
    for j = 1:count % iteration through points        
        sub_magn = magn_map(floor(new_points(j,1) - 3*new_points(j,3)): floor(new_points(j,1) + 3*new_points(j,3)),...
            floor(new_points(j,2) - 3*new_points(j,3)): floor(new_points(j,2) + 3*new_points(j,3)));
        [length,~] = size(sub_magn);
        sub_magn = sub_magn.*gaussianKernel(length,1.5*new_points(j,3)); % gaussian weighting
        sub_orient = orient_map(floor(new_points(j,1) - 3*new_points(j,3)): floor(new_points(j,1) + 3*new_points(j,3)),...
            floor(new_points(j,2) - 3*new_points(j,3)): floor(new_points(j,2) + 3*new_points(j,3)));
        %(subMagn,subOrient,numBin)
        orient_list = generateHist(sub_magn,sub_orient,36);
        [~,numDir] = size(orient_list);
        for k = 1:numDir
            points(end+1,1:5) = new_points(j,1:5);
            points(end,6) = orient_list(k);
        end
    end
    pointLoc = points(:,1:2);
    scales = points(:,3);
    cornerScore = points(:,4);
    scaleScore = points(:,5);
    orientation = points(:,6);
end

function h = gaussianKernel(size,sigma)
    ind = -floor(size/2) : floor(size/2);
    [X Y] = meshgrid(ind, ind);
    h = exp(-(X.^2 + Y.^2) / (2*sigma*sigma));
    h = imresize(h,[size,size]);
    h = h / sum(h(:));
end

function maxOrient = generateHist(subMagn,subOrient,numBin) % eight bins
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
maxOrientList = [];
[maxValue,maxIdx] = max(hist);
maxOrient = binCenter(maxIdx);
if maxOrient >= pi
    maxOrient = maxOrient - 2*pi;
end
maxOrientList = [maxOrientList,maxOrient];
hist(maxIdx) = -1;
[maxValue2,maxIdx2] = max(hist);
if maxValue2 >= 0.8*maxValue
    maxOrientList = [maxOrientList,binCenter(maxIdx2)];
end
hist(maxIdx2) = -1;

[maxValue3,maxIdx3] = max(hist);
if maxValue3 >= 0.8*maxValue
    maxOrientList = [maxOrientList,binCenter(maxIdx3)];
end

maxOrient = maxOrientList;% so it will return a list 
end




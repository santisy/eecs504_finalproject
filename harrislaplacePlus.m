function [pointLoc,scales,cornerScore,scaleScore] = harrislaplacePlus(img,corner_threshold,scale_thresold)
    % modified version dec 08
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
    harris_pts = zeros(0,5);% 4th column: corner score 5th column: laplace response score
   
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
        harris_pts(end+1:end+n,:) = [l,c,repmat(i,[n,1]),max_local(sub2ind([400,400],l,c)),repmat(0,[n,1])];
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
    points = zeros(n,5);
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
        new_points(label,4:5) = points(tdp,4:5);
        new_points(label,3) = min(points(dist_list(:,1)<=2&dist_list(:,2)<=2, 3));
        
        label = label + 1;   
    end
    points = new_points;
    pointLoc = points(:,1:2);
    scales = points(:,3);
    cornerScore = points(:,4);
    scaleScore = points(:,5);

end

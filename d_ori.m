function orient_list = d_ori(points, img)
 s_D = points(3);
 x  = -ceil(3*s_D):ceil(3*s_D);
% if length(x)==1
%     disp('the sigma scale is too small to do Gaussian Derivative')
%     exit();
% end
%dx = x .* exp(-x.*x/(2*s_D*s_D)) ./ (s_D*s_D*s_D*sqrt(2*pi));
%dy = dx';
img = im2double(rgb2gray(img));
dx = [1,0,-1];
dy = dx';

L = imfilter(img,fspecial('gaussian',[ceil(points(3)*6),ceil(points(3)*6)],points(3)),'same');
% image derivatives
Ix = conv2(L, dx, 'same');
Iy = conv2(L, dy, 'same');

magn_map = sqrt(Ix.*Ix + Iy.*Iy);
orient_map = atan2(Iy, Ix);

sub_magn = magn_map(x+points(1),x+points(2));
length_temp = size(sub_magn,1);
sub_magn = sub_magn.*gaussianKernel(length_temp,1.5*s_D); % gaussian weighting
sub_orient = orient_map(x+points(1),x+points(2));
orient_list = generateHist(sub_magn(:),sub_orient(:),36);

end

function h = gaussianKernel(size,sigma)
    ind = -floor(size/2) : floor(size/2);
    [X, Y] = meshgrid(ind, ind);
    h = exp(-(X.^2 + Y.^2) / (2*sigma*sigma));
    h = imresize(h,[size,size]);
    h = h / sum(h(:));
end

function maxOrient = generateHist(subMagn,subOrient,numBin) % eight bins
    % orientation is the orientation ranging from -pi to pi, in rad
    % subMagn is the magnitude
    bins = linspace(-pi,pi,numBin+1);% bin wall
    hist_count = zeros(1,numBin);
    num_of_pix = length(subMagn);
    sub_temp = subOrient*ones(1,numBin+1)- ones(num_of_pix,1)*bins;
    sub_temp(sub_temp<=0) = 0;
    sub_temp(sub_temp>0) = 1;
    index_bin = sum(sub_temp,2);
    for i = 1:numBin
        hist_count(i) = sum(subMagn(index_bin==i));
    end
    [sort_value,sort_index] = sort(hist_count,'descend');
    maxOrient = maxOrientInterp(sort_index(1),hist_count,bins,numBin);
    %Y = circsh
    
    
    
    if sort_value(2)>0.8*sort_value(1)
        maxOrient = [maxOrient maxOrientInterp(sort_index(2),hist_count,bins,numBin)];
    end
end

function maxOrient = maxOrientInterp(sort_index,hist_count,bins,numBin)
    maxOrient = (bins(sort_index)+bins(sort_index+1))/2;
    temp = hist_count;
    temp_index = sort_index;
    shift = 0;
    if (sort_index == 1)
       temp = circshift(hist_count,1,1);
       temp_index = sort_index+1;
       shift = 1;
    end
    if (sort_index == numBin)
       temp = circshift(hist_count,-1,1);
       temp_index = sort_index-1;
       shift = -1;
    end
    xcord = [bins(temp_index-1),bins(temp_index),bins(temp_index+1)];
    y = [temp(temp_index-1),temp(temp_index),temp(temp_index+1)];
    coe = polyfit(xcord,y,2);
    maxOrient = -coe(2)./(2*(coe(1)));
    maxOrient = maxOrient - (2*pi/numBin)*shift; 
end
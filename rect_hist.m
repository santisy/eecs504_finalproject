% extract color histgram in rectangle
function v = rect_hist(im,rect_info,bin)
    % check image range
    if (max(im(:)) >= 1)
        im = single(im)/255;
    end
    [h,w,~] = size(im);
    [x,y] = meshgrid(1:h,1:w);
    x = x(:);
    y = y(:);
    mask = in_rect1([y,x],rect_info);
    mask = reshape(mask,h,w);
    v = histvec(im,mask,bin);
end
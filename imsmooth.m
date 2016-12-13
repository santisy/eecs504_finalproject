function J = imsmooth (I,s)
    [M,N,~] = size(I);
    if (s > 0.01)
        W = ceil(4*s);
        g0 = fspecial('gaussian',[2*W+1,2*W+1],s);
        J = conv2(I,g0,'same');
    else
        J = I;
    end
end
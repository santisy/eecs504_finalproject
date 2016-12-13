function FRAMES = dominantOri(frames, I)
    win_factor = 1.5;
    K = size(frames,2);% get number of keypoints
    FRAMES = [];
    BILINEAR = 0;
    for i = 1:K
        sigma = frames(3,i);
        sigmaw = sigma * win_factor;
        % get blurred image
        im = imsmooth(I,sigma);
        % take derivatives
        dx = 0.5*[-1,0,1]';
        dy = dx';
        Ix = conv2(im,dx,'same');
        Iy = conv2(im,dy,'same');    
        MAG = Ix.*Ix + Iy.*Iy;
        THETA = atan2(Iy,Ix); % -pi~pi
        %=====================================
        % construct orientation histogram     
        %=====================================
        W = max(floor(3*sigmaw),1); % window radius
        x = frames(1,i); % patch center
        y = frames(2,i); %
        xi = floor(x+0.5);
        yi = floor(y+0.5);
        [M,N] = size(im);
        H = zeros(36,1);
        for rx = max(-W, 1-xi):min(W, N-xi)
            for ry = max(-W, 1-yi):min(W, M-yi)
                dx = rx + xi - x;
                dy = ry + yi - y;
                r2 = dx*dx + dy*dy;
                % limit to a circular window
                if (r2 >= W*W + 0.5) 
                    continue; end
                weight = exp(-r2/(2*sigmaw*sigmaw));
                mag = sqrt(MAG(ry+yi,rx+xi));
                theta = fastmod(THETA(ry+yi,rx+xi));

                fbin = 36*theta/(2*pi);
                if (~BILINEAR)
                    H(mod(ceil(fbin+0.5),36)+1) = ...
                        H(mod(ceil(fbin+0.5),36)+1) + mag*weight;
                else
                    % perform bilinear interpolation
                    bin = floor(fbin - 0.5); % -1~35
                    rbin = fbin - bin - 0.5;
                    H(mod(bin + 36,36)+1) = ...
                                H(mod(bin + 36,36)+1)+(1-rbin)*mag*weight;
                    H(mod(bin+1,36)+1) = H(mod(bin+1,36)+1) + weight*mag;
                end
            end
        end
        % smooth histogram
        H = SmoothHistogram(H);
        % Find strongest peak
        maxh = max(H);
        % Scan over the whole histogram, find peaks
        cnt = 0;
        hp = H(36,1);
        tmpFRAMES = [];
        for bin = 1:36
            h0 = H(bin);
            hm = H(mod(bin,36)+1);
            if (h0 > hp && h0 > hm && h0 > 0.8*maxh)
                % set di=0, then take the center of the bins as orientation
                di = -0.5 * (hp-hm) / (hp+hm-2*h0); 
                th = 2*pi*(bin-1+di+0.5)/36;
                new_frame = [frames(1:3,i);th;h0];
                tmpFRAMES = [tmpFRAMES, new_frame];
                cnt = cnt + 1;
            end
        end
        if (cnt >= 2)
            % Take two strongest peak
            [~,idx] = sort(tmpFRAMES(5,:),'ascend');
            FRAMES = [FRAMES, tmpFRAMES(1:4,idx(1:min(4,cnt)))];
        else if (cnt == 1)
                 FRAMES = [FRAMES,tmpFRAMES(1:4)];
             end
        end
    end
end

function th = fastmod (th)
    while(th<0) th = th + 2*pi; end
    while(th>2*pi) th = th - 2*pi; end
end
function h = SmoothHistogram (h)
    for iter = 1:6
        hp = h(36,1);
        for i = 1:36
            hm = h(mod(i,36)+1);
            newh = (hp + h(i) + hm)/3;
            hp = h(i);
            h(i) = newh;
        end
    end
end
function [D,valid_idx] = SIFTdescriptor(FRAMES,I)
    % BY Lowe's Convention
    NBO = 8; % number of histogram bins
    NBP = 4; % number of spatial bins
    magnif = 3;
    wsigma = NBP/2;
    
    % gradient kernel
    dx = 0.5*[-1,0,1];
    dy = dx';
   
    [M,N,~] = size(I);
    K = size(FRAMES,2);
    % create buffer
    D = zeros(NBO*NBP*NBP, K);
    valid_idx = [];
    for i = 1:K
        H = zeros(NBO*NBP*NBP,1);
        % read keypoints
        sigma_i = FRAMES(3,i);
        x = FRAMES(1,i);
        y = FRAMES(2,i);
        xi = floor(x+0.5);
        yi = floor(y+0.5);
        theta_i = fastmod(FRAMES(4,i)); % radian
        cti = cos(theta_i);
        sti = sin(theta_i);
        % blur image
        im = imsmooth(I,sigma_i);
        % calculate window size
        SBP = magnif * sigma_i;
        W = floor( sqrt(2.0) * SBP * (NBP + 1) / 2.0 + 0.5) ;
        
        % Boundary check
        if (xi < 1 || ...
            xi > N || ...
            yi < 1 || ...
            yi > M)
            continue;
        end
        % record valid point indices
        valid_idx = [valid_idx, i];
        
        % calculate gradient and angle
        Ix = conv2(im,dx,'same');
        Iy = conv2(im,dy,'same');
        MAG = sqrt(Ix.*Ix + Iy.*Iy);
        ANG = atan2(Iy,Ix); % -pi ~ pi
        
        for dxi = max(-W,1-xi) : min(W,N-xi)
            for dyi = max(-W,1-yi) : min(W,M-yi)
                mag = MAG(yi+dyi,xi+dxi);
                ang = ANG(yi+dyi,xi+dxi);
                theta = fastmod(-ang+theta_i);
                rx = xi + dxi - x;
                ry = yi + dyi - y;
                % normalize displacement (in bins unit)
                nxi = (cti*rx + sti*ry)/SBP;
                nyi = (-sti*rx + cti*ry)/SBP;
                nti = NBO*theta/(2*pi); 
                % compute Gaussian weight
                weight_i = exp(-(nxi*nxi+nyi*nyi)/(2*wsigma*wsigma));
                
                % Trilinear Interpolation
                binx = floor(nxi - 0.5);
                biny = floor(nyi - 0.5);
                bint = floor(nti);
                
                rbinx = nxi - (binx+0.5);
                rbiny = nyi - (biny+0.5);
                rbint = nti - bint;
                
                % Distibute current sample into 8 adjacent bins
                for dbinx = 0:1
                    for dbiny = 0:1
                        for dbint = 0:1
                            if (binx + dbinx >= -(NBP/2) && ...
                                binx + dbinx < NBP/2 && ...
                                biny + dbiny >= -(NBP/2) && ...
                                biny + dbiny < NBP/2)
                                
                                val = weight_i * mag * ...
                                    abs(1 - dbinx - rbinx) * ...
                                    abs(1 - dbiny - rbiny) * ...
                                    abs(1 - dbint - rbint);
                                % updateHistgram
                                abx = binx + dbinx + 3; 
                                aby = biny + dbiny + 3;
                                %abt = mod(bint + dbint,NBO)+1;
                                if (bint + dbint >= NBO)
                                    abt = bint + dbint -NBO + 1;
                                else
                                    abt = bint + dbint + 1;
                                end
                                histidx = (aby + (abx-1)*NBP - 1) * NBO + abt;
                                H(histidx) = H(histidx) + val;
                            end
                        end
                    end
                end
            end
        end
        % Normalize histogram
        H = H/norm(H);
        % Trucate at 0.2
        for bin = 1:NBO*NBP*NBP
            if (H(bin) > 0.2)
                H(bin) = 0.2;
            end
        end
        % Normalize histogram
        H = H/norm(H);
        % Save
        D(:,i) = H;
    end
    % This function updates Orientation histogram
    % @param rbx 
    %       bin offset along x axis, -2~1
    % @param rby
    %       bin offset along y axis, -2~1
    % @param rbt
    %       bin number of orientaion histogram, 0~8
    % @param updateVal
    %       histogram increment
    % @param hitogram
    %       orientation histogram
    % bin coordinates
    %   (-2,-2), (-1,-2), (0,-2), (1,-2)
    %   (-2,-1), (-1,-1), (0,-1), (1,-1)
    %   (-2, 0), (-1, 0), (0, 0), (1, 0)
    %   (-2, 1), (-1, 1), (0, 1), (1, 1)
    % Inline
    function histogram = updateHistogram(rbx, rby, rbt, updateVal, histogram)
    end   

    function th = fastmod (th)
        while(th<0) th = th + 2*pi; end
        while(th>2*pi) th = th - 2*pi; end
    end
end





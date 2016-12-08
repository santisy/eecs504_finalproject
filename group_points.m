% group points in tha same rect
% return index of points and corresponding rectangle index
function table = group_points(imageid,rect_record,F)
    rectlist = rect_record(rect_record(:,9) == imageid,:);
    rect_num = size(rectlist,1);
    [d,fp_num] = size(F);
    table = [];
    for i = 1:rect_num
        mask = in_rect1([F(2,:);F(1,:)]',repmat(rectlist(i,:),fp_num,1));
        num = sum(mask);
        idx = find(mask == 1);
        % table = [points index, rectangle indx]
        table = [table;idx, i*ones(num,1)];
        % mask F
        F = F.* repmat(~mask',d,1);
    end
end
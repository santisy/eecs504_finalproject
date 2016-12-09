function [matches,distance] = fastmatch(D1,D2,thresh)
    if (nargin == 2)
        thresh = 1;
    end
    costMat = dist2(D1',D2');
    [assignment,~] = munkres(costMat);
    % filter out non distinguishing matching
    match_idx = find(assignment ~= 0)';
    for i = 1:length(match_idx) 
        dist = costMat(match_idx(i),:);
        dist = sort(dist,'ascend');
        if (dist(2)/dist(1) < thresh)
            assignment(match_idx(i)) = 0;
        end
    end
    % finalize matches
    matches(1,:) = find(assignment ~= 0)';
    matches(2,:) = assignment(matches(1,:))';
    idx = sub2ind(size(costMat),matches(1,:)',matches(2,:)');
    distance = costMat(idx);
end
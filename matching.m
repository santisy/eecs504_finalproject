% This script haven't been tested yet!
function [optScore,seg1,seg2] = matching(featureCell1,featureCell2,topN,image1,image2,siftWeight,histWeight,seperationColumn)
% This script visulualize the matching between feature points in each images
% and calculate the most likely correspondance
% expected data input
% featureCell1: a cell structure, each matrix in the cell represents an object(specified by rectangular area)
% Expected structure of each matrix: 1st:2nd columns: column and row coordinate of image point 3:end sift descriptor and histogram discriptor
% topN? number of point-pairs taken into comparison in each matching
% siftWeight: weigthing given to sift
% histWeight: weighting given to histogram
[~,m1] = size(featureCell1);
[~,m2] = size(featureCell2);

optScore = inf
for i = 1:m1
    for j = 1:m2
        [pointsIn1,pointsIn2,currentScoreTemp] = patchMatchScore(featureCell1{i},featureCell2{j},topN,seperationColumn,siftWeight,histWeight);
        if currentScore<=optScoreTemp
            optScoreTemp = currentScore;
            optPoints1 = pointsIn1;
            optPoints2 = pointsIn2;
            seg1Temp = i;
            seg2Temp = j;
        end
    end
end
subplot 221
imshow(image1)
hold on
plot(optPoints1,'r*')
subplot 222
imshow(image2)
hold on
plot(optPoints2,'b*')
seg1 = seg1Temp;
seg2 = seg2Temp;
optScore = optScoreTemp;

function [matchedPointsIn1,matchedPointsIn2,overallScore] = patchMatchScore(featMat1,featMat2,topN,seperationColumn,siftWeight,histWeight)
% This function calculate the matching score between two patches
% topN constrols how many top candicates do we take into consideration
[~,m1] = size(featMat1);
[~,m2] = size(featMat2);
scoreList = zeros(m1,m2)
scoreIdx = zeros(m1,m2,4)
for i = 1:m1
    for j = 1:m2
        sift1 = featMat1(i,3:seperationColumn);
        sift2 = featMat1(j,3:seperationColumn);
        hist1 = featMat1(i,seperationColumn+1:end);
        hist2 = featMat1(i,seperationColumn+1:end);
        score = siftWeight*norm(sift1-sift2)+histWeight*norm(sift1-sift2);
        scoreList(i,j) = score;
        scoreIdx(i,j,1:2) = featMat1(i,1:2);
        scoreIdx(i,j,3:4) = featMat2(j,1:2);
    end
end
matchedPointsIn1 = zeros(topN,2);
matchedPointsIn2 = zeros(ropN,2);
sumScore = [];
for i = 1:topN
    [idx,score] = min(scoreList(:));
    [row,col] = ind2sub(idx);
    matchedPointsIn1(i,:) = scoreIdx(row,col,1:2);
    matchedPointsIn2(i,:) = scoreIdx(row,col,1:2);
    sumScore = [sumScore,scoreList(row,col)];
    scoreList(i,:) = inf;
    scoreList(:,j) = inf;
end
overallScore = sum(sumScore);
    
    



        
        
        






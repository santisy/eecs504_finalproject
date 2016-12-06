# eecs504_finalproject
1. data has detected key points from each card and a circle range.
figs have plots which show detected key points on each card picture. 
2. **harrislaplace.m** gives detected points using multi-scale harris method. return detected points and characteristic
scales. This code is based on [Harris Laplace](https://www.mathworks.com/matlabcentral/fileexchange/17894-keypoint-extraction/content/keypointExtraction/kp_harrislaplace.m).
3. RUN **binary_rect.m** to see the 'rectangle' results. This .m files can generate rectangle ranges, each card picture has 8 rectangles. This file utilizes basic mophological operation, PCA and simple geometry knowledge. The variable _rect_record_ records all rectangles' information, each row represents a rectangle, totally 9 properties:
[x, y, length/2, length_direction, width/2, width_direction,pic#], in which [x,y] is the centre point of a rectangle.
4. **in_rect.m** file can test if certain point is in a given rectangle range.

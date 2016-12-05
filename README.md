# eecs504_finalproject
1. data has detected key points from each card and a circle range.
figs have plots which show detected key points on each card picture. 
2. **harrislaplace.m** gives detected points using multi-scale harris method. return detected points and characteristic
scales.
3. RUN **bianry_rect.m** to see the 'rectangle' results. This .m files can generate rectangle ranges, each card picture has 8 rectangles. The variable _rect_record_ records all rectangles' information, each row represents a rectangle, totally 9 properties:
[x, y, length/2, length_direction, width/2, width_direction,pic#], in which [x,y] is the centre point of a rectangle.
4. **in_rect.m** file can test if certain point is in a given rectangle range.

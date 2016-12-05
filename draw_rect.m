function draw_rect(rect_info)
% about the rectangular
% [x, y, length/2, length_direction, width/2, width_direction,pic#]
syms x y
% here use symbolic mehtod to get the four angle coordinates of each
% rectangle
for i = 1:size(rect_info,1)
    four_corner = zeros(4,2);
    c1 = rect_info(i,3) + rect_info(i,1)*rect_info(i,8) - rect_info(i,2)*rect_info(i,7);
    c2 = c1 - 2*rect_info(i,3);
%     c3 = rect_info(i,6) - rect_info(i,1:2)*rect_info(i,4:5)';
%     c4 = c3 - 2*rect_info(i,6);

    sol1 = solve((x-rect_info(i,1))^2+(y-rect_info(i,2))^2==(rect_info(i,3)^2+rect_info(i,6)^2)...
            ,-rect_info(i,8)*x+rect_info(i,7)*y+c1==0,x,y);
    sol2 = solve((x-rect_info(i,1))^2+(y-rect_info(i,2))^2==(rect_info(i,3)^2+rect_info(i,6)^2)...
            ,-rect_info(i,8)*x+rect_info(i,7)*y+c2==0,x,y);
    four_corner(1:2,1) = double(sol1.x);
    four_corner(1:2,2) = double(sol1.y);
    four_corner(3:4,1) = double(sol2.x);
    four_corner(3:4,2) = double(sol2.y);
%     
    plot(rect_info(:,2), rect_info(:,1),'b+','LineWidth',2)
    plot(four_corner(1:2,2),four_corner(1:2,1),'b-','LineWidth',2);
    plot(four_corner(3:4,2),four_corner(3:4,1),'b-','LineWidth',2);
    plot(four_corner([1,3],2),four_corner([1,3],1),'b-','LineWidth',2);
    plot(four_corner([2,4],2),four_corner([2,4],1),'b-','LineWidth',2);
end
function [nx, ny] = NORMALPROJECTION(BI_x,BI_y,center_x,center_y)

% returns coordinates (x,y) of image point
% find image point x,y by reflection on the boundary

n = sqrt((BI_x-center_x)^2 + (BI_y-center_y)^2);
nx = (BI_x-center_x)/n;
ny = (BI_y-center_y)/n;

return;
end

function [im_x, im_y] = IMAGEPOINT(ghost_x,ghost_y,center_x,center_y,radiusInner,radiusOuter)

% returns coordinates (x,y) of image point
% find image point x,y by reflection on the boundary

ghost_r = sqrt( (ghost_x - center_x)^2 + (ghost_y - center_y)^2 );
ghost_theta = atan2( ghost_y - center_y, ghost_x - center_x );


if (ghost_r <= radiusInner)
    im_r = 2 * radiusInner - ghost_r;
else
    im_r = 2 * radiusOuter - ghost_r;
end


im_theta = ghost_theta;

im_x = center_x + im_r * cos(im_theta);
im_y = center_y + im_r * sin(im_theta);

return;
end

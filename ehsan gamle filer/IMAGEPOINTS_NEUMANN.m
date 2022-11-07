function [im1_x, im1_y, im2_x, im2_y] = IMAGEPOINTS_NEUMANN(ghost_x,ghost_y,center_x,center_y,radiusInner,radiusOuter, h)

% returns coordinates (x,y) of image points

ghost_r = sqrt( (ghost_x - center_x)^2 + (ghost_y - center_y)^2 );
ghost_theta = atan2( ghost_y - center_y, ghost_x - center_x );
im_theta = ghost_theta;

% 
if (ghost_r <= radiusInner)
    im_r = 2 * radiusInner - ghost_r;
else
    im_r = 2 * radiusOuter - ghost_r;
end

im_r = im_r + 0.7 * h;
im1_x = center_x + im_r * cos(im_theta);
im1_y = center_y + im_r * sin(im_theta);

im_r = im_r - 0.7 * h;
im_r = im_r - 0.7 * h;

if (true) % 
end



return;
end

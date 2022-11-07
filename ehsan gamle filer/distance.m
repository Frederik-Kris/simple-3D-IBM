function dl = distance(g_x,g_y,im_x,im_y)

% returns distance between image point and ghost point
% find image point x,y by reflection on the boundary

dl = sqrt((g_x-im_x)^2 + (g_y-im_y)^2);

return;
end

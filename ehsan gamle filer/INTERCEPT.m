function [bi_x, bi_y] = INTERCEPT(x,y,center_x,center_y,radiusInner,radiusOuter)

% returns coordinates (x,y) of intercept point

theta = atan2( y - center_y, x - center_x );
r = sqrt( (x-center_x)^2 + (y-center_y)^2);

if (abs(r-radiusInner)<abs(r-radiusOuter))
    bi_x = center_x + radiusInner * cos(theta);
    bi_y = center_y + radiusInner * sin(theta);
else
    bi_x = center_x + radiusOuter * cos(theta);
    bi_y = center_y + radiusOuter * sin(theta);
end    

return;
end

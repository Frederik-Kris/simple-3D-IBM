function [interp_i, interp_j] = INTERPOLATIONPOINTS(im_x, im_y, x, y, i, j) % i, j: indices of ghost point

global flag;

g_x = x(i,j);
g_y = y(i,j);


if (abs(g_x - im_x) < 1e-12 && abs(g_y - im_y) < 1e-12)  % special case: ghost point is image point

    if (flag(i+1,j+1) < 0)
        interp_i = [ i, i + 1 ];
        interp_j = [ j, j + 1 ];
    elseif (flag(i-1,j-1) < 0)
        interp_i = [ i, i - 1 ];
        interp_j = [ j, j - 1 ];
    else
        error('Unable to create interpolation points.');
    end

else

    % note: assume cartesian grid

    % find i-index
    for i = 1:size(x,1)
        if (x(i,1) >= im_x)
            interp_i_temp = i;
            break;
        end
    end

    % find j-index
    for j = 1:size(x,2)
        if (y(1,j) >= im_y)
            interp_j_temp = j;
            break;
        end
    end

    interp_i = [interp_i_temp - 1, interp_i_temp];
    interp_j = [interp_j_temp - 1, interp_j_temp];
end

return

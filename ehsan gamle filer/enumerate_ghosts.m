function nGhosts = enumerate_ghosts(N,x,y,xcentroidcyl,ycentroidcyl,radiusInner,radiusOuter,c)

global flag;
nGhosts = 0;

% 1. map inside() to flag domain
% 2. check flag directly

for i=1:N
    for j = 1:N
        inside = (x(i,j) - xcentroidcyl)^2 + (y(i,j) - ycentroidcyl)^2 < radiusInner^2 || ...
                 (x(i,j) - xcentroidcyl)^2 + (y(i,j) - ycentroidcyl)^2 > radiusOuter^2 ;
        flag(i,j) = -inside;
    end
end

for i=2:N-1
    for j=2:N-1
        inside_0  = flag(i,j);
        inside_x1 = flag(i+c(1),j);
        inside_x2 = flag(i+c(2),j);
        inside_y1 = flag(i,j+c(1));
        inside_y2 = flag(i,j+c(2));
        if inside_0 && (~inside_x1 || ~inside_x2 || ~inside_y1 || ~inside_y2)
            nGhosts = nGhosts  +  1;
            flag(i,j) = nGhosts; % does not change truth value
        end
    end
end

% swap 0 and -1: fluid means -1. don't touch ghost points
% this is to use flag for equation indexing
for i=1:N
    for j=1:N
        if flag(i,j) == -1
            flag(i,j) = 0;
        elseif flag(i,j) == 0
            flag(i,j) = -1;
        end
    end
end

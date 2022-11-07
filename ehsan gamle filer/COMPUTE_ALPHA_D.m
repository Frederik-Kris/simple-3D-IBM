function [alpha_D, beta_D] = COMPUTE_ALPHA_D(x,y,flag,interpp_i,interpp_j,m,im_x,im_y,...
    b_intercept_x, b_intercept_y)

V_image = [1 im_x im_y im_x*im_y];
Q = ones(1,4);
Vandermonde = zeros(4,4);
for n = 1:4
    
    i = interpp_i(m,n);
    j = interpp_j(m,n);
    
    m_inside = 0;
    if (flag(i,j) > 0)  % (i,j) is a ghost point
        m_inside = flag(i,j);
    end
    
    if m_inside == 0
        int_x = x(interpp_i(m,n), interpp_j(m,n));
        int_y = y(interpp_i(m,n), interpp_j(m,n));
    else
        int_x = b_intercept_x(m_inside);
        int_y = b_intercept_y(m_inside);
        Q(n) = -1;
    end
    
    Vandermonde(n,:) = [1 int_x int_y int_x*int_y];
end

Vand_inv = inv(Vandermonde');

compAlpha = (Vand_inv*V_image');

alpha_D = zeros(1,4);
beta_D = 0;
for n = 1:4
    if Q(n) > 0
        alpha_D(n) = compAlpha(n);
    else
        beta_D = beta_D  +  compAlpha(n);
    end
end

return
end

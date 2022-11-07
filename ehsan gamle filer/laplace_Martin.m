function [l2err, temp, errorDist] = laplace_Martin(N, BCTypeOuter, BCTypeInner, sigma, Tinner, Touter, qNeumann, xCenter, yCenter, r_inner, r_outer, doPlot)

% some defaults so we can call e.g. laplace(100) 

if ~exist('BCTypeOuter','var'); BCTypeOuter = 'Dirichlet'; end
if ~exist('BCTypeInner','var'); BCTypeInner = 'Neumann'  ; end
if ~exist('sigma','var');       sigma       = -1.0       ; end
if ~exist('Tinner','var');      Tinner      = 12         ; end
if ~exist('Touter','var');      Touter      = 6          ; end
if ~exist('qNeumann','var');    qNeumann    = -2.53      ; end
if ~exist('doPlot','var');      doPlot      = true       ; end

% solve laplace equation with dirichlet boundary conditions

global flag gEq;


diamcylInner = 2*r_inner;          % Diameter of cylinder
diamcylOuter = 2*r_outer;          % Diameter of cylinder

% Sanity check

if (strcmp(BCTypeOuter, 'Neumann') && strcmp(BCTypeInner, 'Neumann'))
    error('Neumann-Neumann condition not supported.')
end
if (~strcmp(BCTypeOuter, 'Neumann') && ~strcmp(BCTypeOuter, 'Dirichlet'))
    error('Unknown outer boundary condition type.')
end
if (~strcmp(BCTypeInner, 'Neumann') && ~strcmp(BCTypeInner, 'Dirichlet'))
    error('Unknown inner boundary condition type.')
end


disp('=====================================================');
if (sigma == 0); problemType = 'Laplace'; else problemType = 'Poisson'; end
disp(sprintf('= %s problem, N = %d, %s-%s', problemType, N, BCTypeInner, BCTypeOuter));
disp('=====================================================');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('initialization... '),tic

lx   = 1;              % Length of computational domain in x direction
ly   = 1;              % Length of computational domain in y direction
xcentroidcyl = xCenter;    % x position of centroid cylinder
ycentroidcyl = yCenter;    % y position of centroid cylinder

% Define parameters
radiusInner    = 0.5*diamcylInner;
radiusOuter    = 0.5*diamcylOuter;

% Store mesh
x=zeros(N,N);
y=zeros(N,N);
for j=1:N
    for i=1:N
        x(i,j) = ((i-1)/(N-1))*lx;
        y(i,j) = ((j-1)/(N-1))*ly;
    end
end

h = x(2,1)-x(1,1);

%Compute flag array
%
c = [-1,1];

flag = zeros(N,N);

%figure(1);

toc,fprintf('emumerate ghosts... '),tic

nGhosts = enumerate_ghosts(N,x,y,xcentroidcyl,ycentroidcyl,radiusInner,radiusOuter,c);

ghost_i = zeros(nGhosts,1);
ghost_j = zeros(nGhosts,1);
image_x = zeros(nGhosts,1);
image_y = zeros(nGhosts,1);
delta_l = zeros(nGhosts,1);
n_x = zeros(nGhosts,1);
n_y = zeros(nGhosts,1);
b_intercept_x = zeros(nGhosts,1);
b_intercept_y = zeros(nGhosts,1);
interpp_i = zeros(nGhosts,4);
interpp_j = zeros(nGhosts,4);
int_x = zeros(nGhosts,4);
int_y = zeros(nGhosts,4);
alpha = zeros(nGhosts,4);
beta  = zeros(nGhosts,1);

toc,fprintf('analyze domain... '),tic

for j=1:N
    for i=1:N
       
        if flag(i,j) > 0 % ghost

            m = flag(i,j);

            ghost_i(m) = i;
            ghost_j(m) = j;

            gx = x(i,j);
            gy = y(i,j);
            
            rGhost = sqrt((gx-xcentroidcyl)^2 + (gy-ycentroidcyl)^2);
            if (abs(rGhost - radiusInner) < abs(rGhost - radiusOuter)) % inner boundary
                boundaryType = BCTypeInner;
            else
                boundaryType = BCTypeOuter;
            end
            
            [im_x, im_y] = IMAGEPOINT(x(i,j),y(i,j),xcentroidcyl,ycentroidcyl,radiusInner,radiusOuter);            

            image_x(m) = im_x;
            image_y(m) = im_y;

            dl = distance(x(i,j),y(i,j),im_x,im_y);
            delta_l(m) = dl;

            [bi_x, bi_y] = INTERCEPT(x(i,j),y(i,j),xcentroidcyl,ycentroidcyl,radiusInner,radiusOuter);
            b_intercept_x(m) = bi_x;
            b_intercept_y(m) = bi_y;
            %
            [nx, ny] = NORMALPROJECTION(b_intercept_x(m),b_intercept_y(m),xcentroidcyl,ycentroidcyl);


            if (abs(rGhost - radiusInner) < abs(rGhost - radiusOuter)) % inner boundary
                n_x(m) = nx;
                n_y(m) = ny;
            else % outer boundary
                n_x(m) = -nx;
                n_y(m) = -ny;
            end
            
            [interp_i, interp_j] = INTERPOLATIONPOINTS(im_x, im_y, x, y, i, j);

            interpp_i(m,1) = interp_i(1);
            interpp_i(m,2) = interp_i(1);
            interpp_i(m,3) = interp_i(2);
            interpp_i(m,4) = interp_i(2);
            interpp_j(m,1) = interp_j(1);
            interpp_j(m,2) = interp_j(2);
            interpp_j(m,3) = interp_j(1);
            interpp_j(m,4) = interp_j(2);
                

        end
    end
end

toc,fprintf('compute alpha... '),tic

% loop over ghost points
for m = 1:nGhosts
    % image point
    im_x = image_x(m);
    im_y = image_y(m);

    gx = x( ghost_i(m), ghost_j(m) );
    gy = y( ghost_i(m), ghost_j(m) );
    rGhost = sqrt((gx-xcentroidcyl)^2 + (gy-ycentroidcyl)^2);
    if (abs(rGhost - radiusInner) < abs(rGhost - radiusOuter))
        boundaryType = BCTypeInner;
    else
        boundaryType = BCTypeOuter;
    end
    
    for c = 1:4
        % 4 interpolation points
        int_x(m,c) = x(interpp_i(m,c), interpp_j(m,c));
        int_y(m,c) = y(interpp_i(m,c), interpp_j(m,c));
    end

    if (strcmp(boundaryType, 'Dirichlet'))
        [alpha_, beta_] = COMPUTE_ALPHA_D(x,y,flag,interpp_i,interpp_j,m,im_x,im_y,...
                      b_intercept_x, b_intercept_y);
    else
        [alpha_, beta_] = COMPUTE_ALPHA_N(x,y,flag,interpp_i,interpp_j,m,im_x,im_y,...
            n_x, n_y, b_intercept_x, b_intercept_y);
    end

    alpha(m,:) = alpha_;
    beta(m)  = beta_;
end

toc,fprintf('build matrix... '),tic

% number of nonzero elements in the matrix:
%   5 points for the laplace stencil
%   5 points for the ghost point interpolation
%   minus all ghost points that are part of their own stencil
nEq = sum(sum(flag < 0)) + nGhosts; % fluid points and ghost points
nGhostType23 = sum(sum(abs(alpha) == 0));
NNZ = 5*nEq - nGhostType23; % subtracting all ghost points that are part of their own stencil
II = zeros(NNZ,1);
JJ = zeros(NNZ,1);
SS = zeros(NNZ,1);
temp = 1.0/h^2;
acc_index = 0;
b = zeros(nEq,1);        % right hand side

gEq = zeros(nGhosts,1);  % ghost equation id


% calculate equation indices flag(i,j) < 0: equation for this point is
% the negative of flag(i,j). flag(i,j) > 0: ghost point id. equation index
% is gEq(flag(i,j))

for i = 1:N
    for j = 1:N
        if flag(i,j) < 0
            flag(i,j) = flag(i,j) - acc_index;
            acc_index = acc_index + 1;
        elseif flag(i,j) > 0
            gEq(flag(i,j)) = acc_index + 1;
            acc_index = acc_index + 1;
        end
    end
end

% set matrix elements for interior points

m = 0;
for i = 1:N
    for j = 1:N
        
        if (flag(i,j) >= 0); continue; end

        iEq = -flag(i,j); if (iEq < 0); iEq = gEq(-iEq); end

        b(iEq) = sigma;  % poisson equation            

        m = m + 1;
        II(m) = iEq;
        JJ(m) = iEq;
        SS(m) = -4 * temp;

        m = m + 1;
        jEq = -flag(i,j-1); if (jEq < 0); jEq = gEq(-jEq); end
        II(m) = iEq;
        JJ(m) = jEq;
        SS(m) = temp;

        m = m + 1;
        jEq = -flag(i,j+1); if (jEq < 0); jEq = gEq(-jEq); end
        II(m) = iEq;
        JJ(m) = jEq;
        SS(m) = temp;

        m = m + 1;
        jEq = -flag(i-1,j); if (jEq < 0); jEq = gEq(-jEq); end
        II(m) = iEq;
        JJ(m) = jEq;
        SS(m) = temp;

        m = m + 1;
        jEq = -flag(i+1,j); if (jEq < 0); jEq = gEq(-jEq); end
        II(m) = iEq;
        JJ(m) = jEq;
        SS(m) = temp;
            
    end
end

toc,fprintf('set ghost points... '),tic

for g = 1:nGhosts

    iEq = gEq(g);  % equation for this ghost point
    i = ghost_i(g);
    j = ghost_j(g);
    
    gx = x(i,j);
    gy = y(i,j);
    
    % check closest boundary
    rGhost = sqrt((gx-xcentroidcyl)^2 + (gy-ycentroidcyl)^2);
    if (abs(rGhost - radiusInner) < abs(rGhost - radiusOuter))
        
        boundaryType = BCTypeInner;
        if (strcmp(boundaryType, 'Dirichlet'))
            phi_B = Tinner;
        else
            phi_B = qNeumann;
        end
        
    else
        
        boundaryType = BCTypeOuter;
        if (strcmp(boundaryType, 'Dirichlet'))
            phi_B = Touter;
        else
            phi_B = qNeumann;
        end
        
    end

    if (strcmp(boundaryType, 'Dirichlet'))
        coef = 2 - beta(g);
    else
        coef = delta_l(g) - beta(g);
    end
    
    % central ghost
    m = m + 1;
    II(m) = iEq;
    JJ(m) = iEq;
    SS(m) = 1.0;
    b(iEq) = coef * phi_B;

    for c = 1:4
        
        if (alpha(g,c) == 0)
            continue;
        end

        % interpolation point
        iInt = interpp_i(g,c);
        jInt = interpp_j(g,c);

        jEq = -flag(iInt,jInt); if (jEq < 0); jEq = gEq(-jEq); end
        
        m = m + 1;
        II(m) = iEq;
        JJ(m) = jEq;
        SS(m) = alpha(g,c);

    end
    
end

A = sparse(II,JJ,SS,nEq,nEq,NNZ);


toc,fprintf('solve... '),tic

T = A\b;

toc,fprintf('extract data... '),tic


% analytical solution - depends on BC type
if (strcmp(BCTypeOuter, 'Dirichlet') && strcmp(BCTypeInner, 'Dirichlet'))
    Touter_ = Touter - sigma * radiusOuter^2/4;
    Tinner_ = Tinner - sigma * radiusInner^2/4;
    Bcoef = (log(radiusInner)*Touter_ - log(radiusOuter)*Tinner_) / log(radiusInner/radiusOuter);
    Acoef = 1/log(radiusInner)*(Tinner_ - Bcoef);
else % D-N or N-D
    if (strcmp(BCTypeOuter, 'Dirichlet'))
        rDirichlet = radiusOuter;
        TDirichlet = Touter;
        rNeumann = radiusInner;
        fac = -1;
    else
        rDirichlet = radiusInner;
        TDirichlet = Tinner;
        rNeumann = radiusOuter;
        fac = 1;
    end
    Acoef = rNeumann * fac * qNeumann - sigma * rNeumann^2 / 2;
    Bcoef = TDirichlet - sigma * rDirichlet^2 / 4 - Acoef * log(rDirichlet);
end
% particular solution:  T(r) = sigma * r^2/4 + A ln(r) + B, A and B depend
% on boundary condition type and value
Tanalytic_ = @(r2) sigma*r2/4 + Acoef * log(r2) / 2 + Bcoef;


n_points = 0;
l2err = 0;

temp = zeros(N,N) - NaN;
Tanalytic = zeros(N,N) - NaN;

for i=1:N
    for j=1:N
        if flag(i,j) >= 0; continue; end
        iEq = -flag(i,j);
        temp(i,j) = T(iEq);
        r2 =  (x(i,j) - xcentroidcyl)^2 + (y(i,j) - ycentroidcyl)^2 ;
        Tanalytic(i,j) = Tanalytic_(r2);
        n_points = n_points + 1;
        l2err = l2err + (temp(i,j) - Tanalytic(i,j))^2;
    end
end

l2err = sqrt(l2err/n_points);
errorDist = temp-Tanalytic;
 
if (doPlot)
    xlow  = 0;
    xhigh = lx;
    ylow = 0;
    yhigh = ly;

    close all
    figure(1)
    subplot(2,2,1)
    surface(x,y,temp);
    axis equal;
    axis([xlow xhigh ylow yhigh]);
    colorbar
    shading interp
    title('Numerical solution')

	subplot(2,2,2)
    surface(x,y,Tanalytic-temp);
    axis equal;
    axis([xlow xhigh ylow yhigh]);
    colorbar
    shading interp
    title(sprintf('Absolute error, ||e||_2 = %e', l2err))

	subplot(2,2,3)
    surface(x,y,Tanalytic);
    axis equal;
    axis([xlow xhigh ylow yhigh]);
    colorbar
    shading interp
    title('Analytic solution')
    
    figure(2);
    surface(x,y,temp-Tanalytic);
    axis equal;
    axis([xlow xhigh ylow yhigh]);
    colorbar
    shading interp
    title(sprintf('Absolute error, ||e||_2 = %e', l2err))
    
    figure(3);
    [colors, contourPlot] = contourf(x, y, temp-Tanalytic, 20);
    axis equal;
    %set(contourPlot, 'LineColor', 'none');
    colorbar;
    title(sprintf('Absolute error, ||e||_2 = %e', l2err))
    xlabel('x');
    ylabel('y');
    

end

toc

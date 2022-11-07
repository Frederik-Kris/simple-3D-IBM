
% Update temperature values in the ghost nodes, by first computing the values in the image points:

for n=1:length(ghostNodesIndices) % Loop over ghost/image points
    % Find indices to the 8 surrounding nodes:
    iPrev = floor(xImagePoints(n)/dx)+1;
    jPrev = floor(yImagePoints(n)/dy)+1;
    kPrev = floor(zImagePoints(n)/dz)+1;
    iNext = iPrev + 1;
    jNext = jPrev + 1;
    kNext = kPrev + 1;

    allNeighborsFluid = true; % We set this to false if one or more neighbors are ghosts. If all are fluid we can use the simplified interpolation.
    T_nb = zeros(8,1); % Temperatures in the 8 nodes surrounding the current image point (nb = neighbor)...
    x_nb = zeros(8,1); % and the positions of those 8 points.
    y_nb = zeros(8,1);
    z_nb = zeros(8,1);
    neumannRows(1:8) = false; % This keeps track of which (if any) rows in the Vandermonde matrix must be adapted for neumann conditions.
    normalVectors = zeros(8,3);
    counter = 0; % 1,2,...,8
    for k=[kPrev,kNext] % Loop through the 8 surrounding neighbors.
        for j=[jPrev,jNext]
            for i=[iPrev,iNext]
                counter = counter + 1;
                if isFluidNode(i,j,k) % If it's a fluid node, just add its temperature to the list of interpolation points.
                    T_nb(counter) = intermSolutn(i,j,k);
                    x_nb(counter) = (i-1)*dx;
                    y_nb(counter) = (j-1)*dy;
                    z_nb(counter) = (k-1)*dz;
                else % Ghost node
                    allNeighborsFluid = false;
                    if abs(distanceToCenter(i,j,k) - innerRadius) < abs(distanceToCenter(i,j,k) - outerRadius) % Inner shell is closest
                        if (bcType == "DirDir" || bcType == "DirNeu") % Dirichlet BC on inner shell
                            T_nb(counter) = T_inner;
                        else
                            neumannRows(counter) = true;
                            T_nb(counter) = tempGradInner;
                        end
                    else % Outer shell is closest
                        if (bcType == "DirDir" || bcType == "NeuDir") % Dirichlet BC on outer shell
                            T_nb(counter) = T_outer;
                        else
                            neumannRows(counter) = true;
                            T_nb(counter) = tempGradOuter;
                        end
                    end
                    if (k==NK && caseType == "Cyl")    % These are special cases, where the surrounding node is neither fluid nor ghost.
                        kAdjust = NK-1;                % We employ this quick fix, which only works for the cylinder case, with periodic
                    elseif (k==1 && caseType == "Cyl") % BC in z-direction.
                        kAdjust = 2;
                    else
                        kAdjust = k;
                    end
                    thisGP_index = ghostIndexMap( get1D_index(i,j,kAdjust,NI,NJ) );       % What's happening here is we first get the 1D index to the
                    x_nb(counter) = xBodyIntercept(thisGP_index);                  % ghost node i,j,k. Then we use the map to find that index's
                    y_nb(counter) = yBodyIntercept(thisGP_index);                  % position in the 'ghostNodesIndices' array. This NOT necessarily
                    z_nb(counter) = zBodyIntercept(thisGP_index) + (k-kAdjust)*dz; % equal to 'n', since multiple ghost nodes can surround the IP.
                    if neumannRows(counter) && caseType=="Cyl"
                        normalVectors(counter,:) = [(i-1)*dx-xCenter, (j-1)*dy-yCenter, 0]; % normal vector if cylinder case
                    elseif neumannRows(counter) && caseType=="Sphere"
                        normalVectors(counter,:) = [(i-1)*dx-xCenter, (j-1)*dy-yCenter, (k-1)*dz-zCenter]; % normal vector if sphere case
                    end
                end
            end
        end
    end
    % At this stage we have the temperature and position of 8 points surrounding the IP.
    % Then we can find the temperature in the IP by simplified 3D interpolation or full trilinear interpolation.
    if allNeighborsFluid
        xPrev = (iPrev-1)*dx; % Coordinates to the point with lowest indices
        yPrev = (jPrev-1)*dy;
        zPrev = (kPrev-1)*dz;
        dx_IP = xImagePoints(n)-xPrev; % ...and distance from that point to the IP
        dy_IP = yImagePoints(n)-yPrev;
        dz_IP = zImagePoints(n)-zPrev;
        T12 = T_nb(1) + dx_IP * (T_nb(2)-T_nb(1)) / dx; % Linear interpolation between point 1 and 2
        T34 = T_nb(3) + dx_IP * (T_nb(4)-T_nb(3)) / dx; % Linear interpolation between point 3 and 4
        T56 = T_nb(5) + dx_IP * (T_nb(6)-T_nb(5)) / dx; % Linear interpolation between point 5 and 6
        T78 = T_nb(7) + dx_IP * (T_nb(8)-T_nb(7)) / dx; % Linear interpolation between point 7 and 8
        T1234 = T12 + dy_IP * (T34-T12) / dy; % ...and then we keep interpolating...
        T5678 = T56 + dy_IP * (T78-T56) / dy;
        T_IP = T1234 + dz_IP * (T5678-T1234) / dz; % ...until we only have one temperature.
    else
        coefficientMatrix = zeros(8,8); % Vandermonde matrix. We make 8 rows for it:
        for m=1:8
            if ~neumannRows(m) % 'Regular' row
                coefficientMatrix(m,:) = [1, x_nb(m), y_nb(m), z_nb(m), x_nb(m)*y_nb(m), x_nb(m)*z_nb(m), y_nb(m)*z_nb(m), x_nb(m)*y_nb(m)*z_nb(m)];
            else % Special 'Neumann' row. First find unit normal. nVec = normalVector
                nVec = normalVectors(m,:);
                if norm(nVec) < innerRadius % Make it point from fluid to solid
                    nVec = -nVec;
                end
                nVec = nVec / norm(nVec); % Make it UNIT normal
                coefficientMatrix(m,:) = [0, nVec(1), nVec(2), nVec(3), x_nb(m)*nVec(2)+y_nb(m)*nVec(1), x_nb(m)*nVec(3)+z_nb(m)*nVec(1), y_nb(m)*nVec(3)+z_nb(m)*nVec(2), x_nb(m)*y_nb(m)*nVec(3)+x_nb(m)*z_nb(m)*nVec(2)+y_nb(m)*z_nb(m)*nVec(1)];
            end
        end
        a = coefficientMatrix \ T_nb; % Solving the linear system gives the coefficients, that we store in 'a'.
        % Then use the coefficients to compute the temperature in the IP as a function of x,y,z:
        T_IP = a(1) + a(2)*xImagePoints(n) + a(3)*yImagePoints(n) + a(4)*zImagePoints(n) + a(5)*xImagePoints(n)*yImagePoints(n) + a(6)*xImagePoints(n)*zImagePoints(n) + a(7)*yImagePoints(n)*zImagePoints(n) + a(8)*xImagePoints(n)*yImagePoints(n)*zImagePoints(n);
    end
    % Use the BC to get the temperature in the actual ghost node:
    rGhost = distanceToCenter(ghostNodesIndices(n)); % Distance from cylinder/sphere center to ghost node.
    if abs(rGhost - innerRadius) < abs(rGhost - outerRadius) % Inner shell is closest
        if (bcType == "DirDir" || bcType == "DirNeu") % Dirichlet on this shell
            intermSolutn(ghostNodesIndices(n)) = 2*T_inner - T_IP;
        else % Neumann on this shell
            normalProbeLength = sqrt((xImagePoints(n)-xBodyIntercept(n))^2+(yImagePoints(n)-yBodyIntercept(n))^2+(zImagePoints(n)-zBodyIntercept(n))^2);
            intermSolutn(ghostNodesIndices(n)) = T_IP + 2*normalProbeLength*tempGradInner;
        end
    else % Outer shell is closest
        if (bcType == "DirDir" || bcType == "NeuDir") % Dirichlet on this shell
            intermSolutn(ghostNodesIndices(n)) = 2*T_outer - T_IP;
        else % Neumann on this shell
            normalProbeLength = sqrt((xImagePoints(n)-xBodyIntercept(n))^2+(yImagePoints(n)-yBodyIntercept(n))^2+(zImagePoints(n)-zBodyIntercept(n))^2);
            intermSolutn(ghostNodesIndices(n)) = T_IP + 2*normalProbeLength*tempGradOuter;
        end
    end
end











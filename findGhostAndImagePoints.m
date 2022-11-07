
% Find the ghost nodes among the solid nodes, and find their respective image points:
counter = 0;
ghostNodesIndices = [];
for n=find(~isFluidNode)' % loop through solid nodes
    [i,j,k] = get3D_indices(n,NI,NJ);  % get the 3D indices to find neighbors
    if i<NI && isActiveNode(i+1,j,k) ...
    || i>1  && isActiveNode(i-1,j,k) ... % if one or more neighbor is an active node...
    || j<NJ && isActiveNode(i,j+1,k) ... % (The shortcutting ANDs here ensure that we
    || j>1  && isActiveNode(i,j-1,k) ... %  don't check out of bounds)
    || k<NK && isActiveNode(i,j,k+1) ...
    || k>1  && isActiveNode(i,j,k-1)
        counter = counter + 1;
        ghostNodesIndices(counter) = n; % ...then we know it's a ghost node.
        ghostNodePosition = [(i-1)*dx, (j-1)*dy, (k-1)*dz]; % The absolute position vector of the ghost node
        % Find the normal probe from ghost node to relevant surface:
        if (caseType == "Cyl")
            radialPosition = [ghostNodePosition(1)-xCenter, ghostNodePosition(2)-yCenter, 0]; % vector from cylinder centroid to ghost node.
        elseif (caseType == "Sphere")
            radialPosition = [ghostNodePosition(1)-xCenter, ghostNodePosition(2)-yCenter, ghostNodePosition(3)-zCenter]; % vector from sphere center to ghost node.
        end
        if distanceToCenter(n) < innerRadius
            normalProbe = radialPosition * (innerRadius - distanceToCenter(n)) / distanceToCenter(n);
        elseif distanceToCenter(n) > outerRadius
            normalProbe = radialPosition * (outerRadius - distanceToCenter(n)) / distanceToCenter(n);
        end
        bodyInterceptPoint = ghostNodePosition + normalProbe;
        xBodyIntercept(counter) = bodyInterceptPoint(1);
        yBodyIntercept(counter) = bodyInterceptPoint(2);
        zBodyIntercept(counter) = bodyInterceptPoint(3);
        imagePointPosition = bodyInterceptPoint + normalProbe;
        xImagePoints(counter) = imagePointPosition(1);
        yImagePoints(counter) = imagePointPosition(2);
        zImagePoints(counter) = imagePointPosition(3);
    end
end
ghostIndexMap = containers.Map(ghostNodesIndices, 1:counter); % Here we create a map so we can find the index to 'ghostNodesIndices'. It becomes like the inverse of looking up an element in the array.

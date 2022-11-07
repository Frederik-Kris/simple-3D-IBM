
% Flagging and finding the indices to different node categories:

distanceToCenter = zeros(NI,NJ,NK); % Find distance between nodes and cylinder centroid or sphere center.
for i=1:NI
    for j=1:NJ
        x=(i-1)*dx;
        y=(j-1)*dy;
        if (caseType == "Cyl")
            distanceToCenter(i,j,:) = sqrt( (x-xCenter)^2 + (y-yCenter)^2 );
        elseif (caseType == "Sphere")
            for k=1:NK
                z=(k-1)*dz;
                distanceToCenter(i,j,k) = sqrt( (x-xCenter)^2 + (y-yCenter)^2 + (z-zCenter)^2 );
            end
        end
    end
end
isFluidNode = distanceToCenter >= innerRadius - epsilon ... % Flag array. True for fluid nodes,
            & distanceToCenter <= outerRadius + epsilon;    % and false for solid nodes.
fluidNodesIndices = find(isFluidNode)';
isSurfaceNode = distanceToCenter >= innerRadius - epsilon ... % Flag array. True for nodes that lie on
              & distanceToCenter <= innerRadius + epsilon ... % the immersed boundary (within tolerance),
              & bcType ~= "NeuDir"                        ... % if that IB has Dirichlet condition.
              | distanceToCenter >= outerRadius - epsilon ... % They are defined as fluid nodes, but get
              & distanceToCenter <= outerRadius + epsilon ... % their values directly from BC at the IB.
              & bcType ~= "DirNeu";
surfaceNodesIndices = find(isSurfaceNode)';
isInteriorNode = zeros(NI,NJ,NK);
isInteriorNode(:,:,:) = false;         % Flag array. True for nodes that are not on the domain boundaries.
isInteriorNode(2:NI-1,2:NJ-1,2:NK-1) = true;
isActiveNode = isInteriorNode & isFluidNode & ~isSurfaceNode; % In active nodes the governing eq is solved
activeNodesIndices = find(isActiveNode)';  % Array of 1D indices to the interior fluid nodes (excluding the boundaries).
if (caseType == "Cyl") % For the cylinder case, find the nodes at boundaries that the cylinder intersects.
    isLowerBoundary = zeros(NI,NJ,NK);
    isLowerBoundary(:,:,:) = false;
    isLowerBoundary(:,:,1) = true;
    lowerIndicesBC = find(isFluidNode & isLowerBoundary)';
    isUpperBoundary = zeros(NI,NJ,NK);
    isUpperBoundary(:,:,:) = false;
    isUpperBoundary(:,:,NK) = true;
    upperIndicesBC = find(isFluidNode & isUpperBoundary)';
end

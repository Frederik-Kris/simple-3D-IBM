
% Show pos of ghost and image points:
function showGhostAndImagePositions(figNum)
figure(figNum);
for n = 1:length(ghostNodesIndices)
    [i,j,k] = get3D_indices(ghostNodesIndices(n),NI,NJ);
    if k==2
        xVals = [(i-1)*dx, xBodyIntercept(n), xImagePoints(n)];
        yVals = [(j-1)*dy, yBodyIntercept(n), yImagePoints(n)];
        plot(xVals,yVals,'ro-');
    end
end
end

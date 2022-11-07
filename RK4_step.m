function slope = RK4_step(T, alpha, activeNodesIndices, NI, NJ, dx, dy, dz)
%RK4_STEP Computes the time slope or residual function of the solution, only in active fluid nodes.
%   Detailed explanation goes here
[i,j,k] = get3D_indices(activeNodesIndices,NI,NJ);
indexM = activeNodesIndices;                % Mid, (i,j,k)
indexE = get1D_index(i+1, j, k, NI, NJ);    % East, i+1
indexW = get1D_index(i-1, j, k, NI, NJ);    % West, i-1
indexN = get1D_index(i, j+1, k, NI, NJ);    % North, j+1
indexS = get1D_index(i, j-1, k, NI, NJ);    % South, j-1
indexU = get1D_index(i, j, k+1, NI, NJ);    % Up, k+1
indexD = get1D_index(i, j, k-1, NI, NJ);    % Down, k-1
slope = alpha * ( (T(indexE)-2*T(indexM)+T(indexW))/dx^2 ...
                + (T(indexN)-2*T(indexM)+T(indexS))/dy^2 ...
                + (T(indexU)-2*T(indexM)+T(indexD))/dz^2 );
end


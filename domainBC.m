% Apply the boundary conditions for the fluid nodes at the boundary of the domain. 
% NOT the immersed boundary.

if caseType == "Cyl"
    [i,j,k] = get3D_indices(lowerIndicesBC,NI,NJ);  % Get 3D indices for fluid nodes at the lower boundary (k=1, z=0)
    upperBoundaryAdjacentIndices = get1D_index(i,j,NK-1, NI,NJ); % Get 1D indices to all fluid nodes adjacent to the upper boundary, k=NK-1
    intermSolutn(lowerIndicesBC) = intermSolutn(upperBoundaryAdjacentIndices); % Periodic BC, take values from boundary adjacent node at upper boundary.
    [i,j,k] = get3D_indices(upperIndicesBC,NI,NJ);  % Get 3D indices for fluid nodes at the upper boundary (k=NK, z=L_z)
    lowerBoundaryAdjacentIndices = get1D_index(i,j,2, NI,NJ); % Get 1D indices to all fluid nodes adjacent to the lower boundary, k=2
    intermSolutn(upperIndicesBC) = intermSolutn(lowerBoundaryAdjacentIndices); % Periodic BC, take values from boundary adjacent node at lower boundary.
end

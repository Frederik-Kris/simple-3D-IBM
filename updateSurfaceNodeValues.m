
% Update values in the surface nodes, the ones (if any) that coincide with an IB:

for n=1:length(surfaceNodesIndices) % Loop over surface nodes
    rNode = distanceToCenter(surfaceNodesIndices(n)); % Distance from cylinder/sphere center to surface node.
    if abs(rNode-innerRadius) <= epsilon
        if (bcType == "DirDir" || bcType == "DirNeu") % Dirichlet on this (inner) IB
            intermSolutn(surfaceNodesIndices(n)) = T_inner;
        end
    elseif abs(rNode-outerRadius) <= epsilon
        if (bcType == "DirDir" || bcType == "NeuDir") % Dirichlet on this (outer) IB
            intermSolutn(surfaceNodesIndices(n)) = T_outer;
        end
    else
        disp("Impossible situation when evaluating surface nodes.");
        fprintf("rNode-innerRadius=%e",rNode-innerRadius);
        fprintf("rNode-outerRadius=%e",rNode-outerRadius);
        fprintf("epsilon=%e",epsilon);
    end
end

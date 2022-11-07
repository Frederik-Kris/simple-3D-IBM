
% Compute analytical solution and store it in an array, so it can be compared to the numerical solution.

% The analytical solution always have two constants A and B, depending on the case:
if caseType == "Cyl"
    if bcType == "DirDir"
        B = (T_outer*log(innerRadius) - T_inner*log(outerRadius)) / log(innerRadius/outerRadius);
        A = (T_inner - B) / log(innerRadius);
    elseif bcType == "DirNeu"
        A = outerRadius * tempGradOuter;
        B = T_inner - A*log(innerRadius);
    elseif bcType == "NeuDir"
        A = innerRadius * -tempGradInner;
        B = T_outer - A*log(outerRadius);
    end
elseif caseType == "Sphere"
    if bcType == "DirDir"
        B = (outerRadius*T_outer - innerRadius*T_inner) / (outerRadius - innerRadius);
        A = innerRadius * (T_inner - B);
    elseif bcType == "DirNeu"
        A = -tempGradOuter * outerRadius^2;
        B = T_inner - A / innerRadius;
    elseif bcType == "NeuDir"
        A = tempGradInner * innerRadius^2;
        B = T_outer - A / outerRadius;
    end
end

% Next, find T_analytical. Set NaN on solid nodes except ghosts:
T_a = NaN * zeros(NI,NJ,NK);
r = distanceToCenter([fluidNodesIndices, ghostNodesIndices]);
if caseType == "Cyl"
    T_a([fluidNodesIndices, ghostNodesIndices]) = A * log(r) + B; % Analytical solution in fluid and ghost nodes, cylinder case.
elseif caseType == "Sphere"
    T_a([fluidNodesIndices, ghostNodesIndices]) = A ./ r + B; % Analytical solution in fluid and ghost nodes, sphere case.
end
errorAbs = T-T_a; % Absolute error distribution
if caseType == "Cyl"
    error1Norm = dx*dy*sum(abs(errorAbs([activeNodesIndices, surfaceNodesIndices])));
    error2Norm = sqrt(dx*dy*sum(errorAbs([activeNodesIndices, surfaceNodesIndices]).^2));
elseif caseType == "Sphere"
    error1Norm = dx*dy*dz*sum(abs(errorAbs([activeNodesIndices, surfaceNodesIndices])));
    error2Norm = sqrt(dx*dy*dz*sum(errorAbs([activeNodesIndices, surfaceNodesIndices]).^2));
end
if testAllCases
    errors_1Norm(currentCaseNumber, currentGridNumber) = error1Norm;
    errors_2Norm(currentCaseNumber, currentGridNumber) = error2Norm;
end
fprintf("2-norm of error: %.3e \n\n",error2Norm);



